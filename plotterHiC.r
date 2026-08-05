#!/usr/bin/env Rscript
# =============================================================================
# plotterHiC.r
#
# Rotated 45-degree Hi-C heatmaps stacked over a single genome-browser panel
# that can carry bigwig signal tracks, TAD/loop arcs and any BED annotation
# (full gene models with exons/introns/chevrons, or plain interval boxes).
#
# Author: Oscar Amaury Aguilar Lomas
# Date:   May 5, 2026
#
# -----------------------------------------------------------------------------
# USAGE
#   Rscript plotterHiC.r config/plotterHiC_mESC_ysEC.yaml
#       -> one figure per entry in global.loci
#
#   Rscript plotterHiC.r config/plotterHiC_mESC_ysEC.yaml chr2:112225000-114225000 Actc1
#       -> single figure; the CLI region and name override global.loci
#
# -----------------------------------------------------------------------------
# DESIGN RULES (why the config looks the way it does)
#
#   1. EVERY path is a full, independent path. Nothing is assembled from a
#      shared directory + a file name, so any input can live anywhere.
#   2. EVERY collection is a NAMED MAP (name: {settings}), never a set of
#      parallel lists that must stay index-aligned.
#   3. EVERY section has a `defaults:` block inherited by its entries; any key
#      set inside an entry overrides it.
#   4. Nothing is hardcoded to a fixed number of items: N Hi-C files, N panels,
#      N tracks of any type.
#   5. ORDER IS DECLARATION ORDER. Panels are drawn top to bottom in the order
#      they appear under `hic.panels`, and tracks in the order they appear
#      under `tracks.list`. There are no y coordinates to keep in sync: you
#      only ever set a track's `height` and the `spacing` between tracks, and
#      the vertical layout is computed from those.
#
# -----------------------------------------------------------------------------
# LAYOUT
#
#   One heatmap panel per `hic.panels` entry, then ONE track panel at the
#   bottom holding every entry of `tracks.list`. Relative panel heights come
#   from each panel's own `height`, and the track panel's from
#   `tracks.panel_height`.
#
#   +--------------------------------------------------+
#   | heatmap  hic.panels[[1]]            height: 1     |
#   | heatmap  hic.panels[[2]]            height: 1     |
#   | heatmap  hic.panels[[3]]  (compare) height: 1     |
#   +--------------------------------------------------+
#   | track panel                tracks.panel_height    |
#   |   ~~~~~   tracks.list[[1]]  type: structure       |  <- first = TOP
#   |   /\_/\   tracks.list[[2]]  type: mark            |
#   |   ->--[]  tracks.list[[3]]  type: annotation      |  <- last  = BOTTOM
#   +--------------------------------------------------+
#
#   Inside the track panel each track owns a horizontal BAND of height
#   `height`, separated from its neighbours by `spacing`. Everything a track
#   draws is expressed relative to its own band, so changing one track's
#   height or reordering the list never requires touching another track.
#
# -----------------------------------------------------------------------------
# HI-C MODEL
#
#   `hic.files` is a NAMED REGISTRY of .hic files. Panels reference those
#   names, so any file can be plotted on its own or compared against any other:
#
#     files:  { ysEC: /path/a.hic, mESC: /path/b.hic, KO: /path/c.hic }
#     panels:
#       - { name: ysEC,          hic: ysEC,           matrix: log2 }
#       - { name: "mESC - ysEC", hic: mESC, vs: ysEC, matrix: subtraction }
#       - { name: "KO / ysEC",   hic: KO,   vs: ysEC, matrix: ratio }
#
#   matrix types:
#     counts       raw contacts divided by the chromosome-wide mean  (needs hic)
#     log2         log2 of the above                                 (needs hic)
#     subtraction  (hic - vs), rescaled by the mean of both means    (needs vs)
#     ratio        hic / vs                                          (needs vs)
#
#   Counts are normalised per file against that file's own chromosome-wide mean
#   contact density, which is what makes libraries of different depth (e.g. a
#   381M vs a 339M read pair experiment) directly comparable.
# =============================================================================

start.time <- Sys.time()
options(scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(patchwork)
  library(GenomicRanges)
  library(rtracklayer)
  library(ggrepel)
  library(yaml)
})
# straw is called qualified (strawr::straw). Install with:
#   remotes::install_github("aidenlab/straw/R")

# =============================================================================
# HELPERS
# =============================================================================

# NULL-coalesce. Also treats a zero-length value as absent, so an empty YAML
# key (`foo:` with nothing after it) falls through to the default.
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# Resolve one option with the precedence: entry > type defaults > fallback.
opt <- function(entry, defaults, key, fallback = NULL) {
  entry[[key]] %||% defaults[[key]] %||% fallback
}

# Same, for a nested block (e.g. `labels:`): merges the entry's sub-list on top
# of the defaults' sub-list, so a track can override just one label setting.
opt_block <- function(entry, defaults, key) {
  modifyList(defaults[[key]] %||% list(), entry[[key]] %||% list())
}

# Effective settings for one entry, resolved in three layers:
#
#   type defaults  <  preset  <  the entry itself
#
# The middle layer exists because one `defaults` block per TYPE is too coarse
# once you keep several tracks of the same subtype: four loop tracks (WT, KO,
# Up, Down) are all `type: structure`, yet they share a look that the TAD
# tracks do not. A named preset is that shared look, written once.
# modifyList recurses, so nested blocks such as `labels:` merge key by key
# instead of one layer replacing the whole block.
merge_settings <- function(entry, type_defaults, presets, what = "entry") {
  base <- type_defaults %||% list()
  pn <- entry$preset
  if (!is.null(pn)) {
    p <- presets[[as.character(pn)]]
    if (is.null(p)) {
      stop(what, ": unknown preset '", pn, "'. Defined presets: ",
           if (length(presets)) paste(names(presets), collapse = ", ") else "(none)",
           call. = FALSE)
    }
    base <- modifyList(base, p)
  }
  modifyList(base, entry)
}

# Read the YAML as raw bytes. `readLines` truncates the file at the first
# non-ASCII byte when the execution locale is C (the default for Rscript
# launched with LANG unset), and yaml::read_yaml then silently returns only the
# keys above that byte. Going through readBin avoids the locale entirely.
read_yaml_utf8 <- function(path) {
  if (!file.exists(path)) stop("config file not found: ", path, call. = FALSE)
  txt <- rawToChar(readBin(path, "raw", file.info(path)$size))
  Encoding(txt) <- "UTF-8"
  yaml::yaml.load(txt)
}

# Accepts chr2:112225000-114225000 and chr2:112225000:114225000, with or
# without the "chr" prefix, and tolerating thousands separators.
parse_region <- function(x) {
  x <- trimws(as.character(x))
  m <- regmatches(x, regexec("^(chr)?([^:]+):([0-9,._]+)[[:space:]]*[-:][[:space:]]*([0-9,._]+)$", x))[[1]]
  if (!length(m)) {
    stop("cannot parse region '", x, "'. Expected chr2:112225000-114225000", call. = FALSE)
  }
  num <- function(s) as.numeric(gsub("[,_ ]", "", s))
  out <- list(chr = m[3], start = num(m[4]), end = num(m[5]))
  if (out$end <= out$start) stop("region '", x, "' ends at or before it starts", call. = FALSE)
  out
}

# A true/false switch, absent meaning TRUE. Anything that is not a real YAML
# boolean is an error rather than a silent FALSE: a typo like `enabled: tru`
# arrives as the string "tru", and quietly dropping the track it was meant to
# turn ON is indistinguishable from the file simply having no features.
as_flag <- function(x, what) {
  if (is.null(x) || length(x) == 0) return(TRUE)
  if (is.logical(x) && length(x) == 1 && !is.na(x)) return(x)
  stop(what, ": 'enabled' must be true or false, got '", paste(x, collapse = " "), "'",
       call. = FALSE)
}

# Numeric option that may legitimately be absent. YAML has no R `NA`, so the
# usual ways of writing "let the data decide" -- NA, null, auto, blank -- all
# arrive as strings; convert them without the coercion warning.
num_or_na <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  if (is.character(x) && tolower(trimws(x)) %in% c("na", "null", "auto", "none", "")) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(x))
}

# A colour-scale limit, written in any of three ways:
#
#   NA / null / auto   the data's own min (or max) in the visible window
#   a number           that exact value, e.g. 3 or -40
#   qN, N in [0,1]     the Nth quantile of the data, e.g. q0.25 = first
#                      quartile, q0.99 = clip the top 1% of outliers
#
# Quantiles matter here because a single saturated bin can otherwise flatten a
# whole map: q0.01/q0.99 gives a robust range without hand-tuning per locus.
parse_limit <- function(x, what) {
  if (is.null(x) || length(x) == 0) return(list(kind = "auto"))
  s <- trimws(as.character(x)[1])
  if (tolower(s) %in% c("na", "null", "auto", "none", "")) return(list(kind = "auto"))
  if (grepl("^[qQ]", s)) {
    p <- suppressWarnings(as.numeric(sub("^[qQ]", "", s)))
    if (is.na(p) || p < 0 || p > 1) {
      stop(what, ": '", s, "' is not a valid quantile. Use qN with N between 0 and 1, ",
           "e.g. q0.25 or q0.99", call. = FALSE)
    }
    return(list(kind = "quantile", value = p))
  }
  v <- suppressWarnings(as.numeric(s))
  if (is.na(v)) {
    stop(what, ": '", s, "' is neither a number, a quantile (qN) nor NA", call. = FALSE)
  }
  list(kind = "fixed", value = v)
}

# Turn a parsed limit into a number, given the values actually being drawn.
resolve_limit <- function(spec, vals, fallback) {
  switch(spec$kind,
    fixed    = spec$value,
    quantile = unname(quantile(vals, spec$value, na.rm = TRUE)),
    fallback(vals, na.rm = TRUE))          # "auto"
}

need_file <- function(path, what) {
  if (is.null(path) || !nzchar(path)) stop(what, ": no file given", call. = FALSE)
  if (!file.exists(path)) stop(what, ": file not found -> ", path, call. = FALSE)
  path
}

# "R,G,B" -> "#RRGGBB". Hex strings and R color names pass through untouched.
rgb_from_string <- function(s) {
  s <- as.character(s)[1]
  if (is.na(s) || !grepl(",", s, fixed = TRUE)) return(s)
  v <- as.integer(strsplit(s, ",", fixed = TRUE)[[1]])
  grDevices::rgb(v[1], v[2], v[3], maxColorValue = 255)
}

# =============================================================================
# GENOMES
# =============================================================================
mm10 <- data.table(
  chr  = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X","Y"),
  size = c(195471971,182113224,160039680,156508116,151834684,149736546,145441459,129401213,124595110,
           130694993,122082543,120129022,120421639,124902244,104043685,98207768,94987271,90702639,
           61431566,171031299,91744698))
mm10_test <- data.table(chr = c("18","19"), size = c(90702639, 61431566))
hg38 <- data.table(
  chr  = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),
  size = c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,
           133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,
           58617616,64444167,46709983,50818468,156040895,57227415))
hg19 <- data.table(
  chr  = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),
  size = c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,
           135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,
           59128983,63025520,48129895,51304566,155270560,59373566))
genome.list <- list(mm10 = mm10, hg38 = hg38, hg19 = hg19, mm10_test = mm10_test)

# =============================================================================
# CONFIG PARSING
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else "config/plotter_config.yaml"
cfg <- read_yaml_utf8(cfg_path)

# Old schemas are refused rather than half-honoured: silently ignoring a
# renamed section would produce a plausible but wrong figure.
# NOTE: `[[` with exact=TRUE, because `$` partial-matches ("file" -> "files").
has_key <- function(x, k) !is.null(x[[k, exact = TRUE]])
if (has_key(cfg, "genes") || has_key(cfg, "marks") || has_key(cfg, "annotations")) {
  stop("this config uses the old 'genes:' / 'marks:' / 'annotations:' sections.\n",
       "  Every track now lives in ONE ordered list, 'tracks.list', where each entry\n",
       "  declares 'type: mark | structure | annotation' and the declaration order is\n",
       "  the order on the figure. See config/plotterHiC_mESC_ysEC.yaml", call. = FALSE)
}
if (!is.null(cfg$hic) && (has_key(cfg$hic, "file") || has_key(cfg$hic, "ctrl_file"))) {
  stop("this config uses the old 'hic.file' / 'hic.ctrl_file' pair. Hi-C files are now a\n",
       "  named registry under 'hic.files:', referenced by 'hic:' / 'vs:' in each panel.\n",
       "  See config/plotterHiC_mESC_ysEC.yaml", call. = FALSE)
}

g <- cfg$global %||% list()

genome_name      <- g$genome %||% "mm10"
resol            <- as.numeric(g$resolution %||% 10000)
padding_fraction <- as.numeric(g$padding_fraction %||% 0.333)

chr.list <- genome.list[[genome_name]]
if (is.null(chr.list)) {
  stop("unknown genome '", genome_name, "'. Available: ",
       paste(names(genome.list), collapse = ", "), call. = FALSE)
}

out_cfg     <- g$output %||% list()
dir_plots   <- out_cfg$dir %||% "."
plot_format <- tolower(out_cfg$formato %||% "png")
plot_width  <- as.numeric(out_cfg$plot_width  %||% 10)   # cm
plot_height <- as.numeric(out_cfg$plot_height %||% 15)   # cm
plot_res    <- as.numeric(out_cfg$dpi %||% 500)          # png only
if (!plot_format %in% c("png", "pdf")) {
  stop("global.output.formato must be 'png' or 'pdf', got '", plot_format, "'", call. = FALSE)
}

# ---- loci: a CLI region wins over global.loci -------------------------------
# A locus entry is EITHER a bare region string OR a braced block carrying the
# region plus per-locus overrides:
#
#   loci: {
#     Compartments: chr2:90000000-120000000,                       # inherits all
#     TADs: { region: chr2:111000000-115000000, resolution: 10000, plot_width: 14 }
#   }
#
# A 30 Mb compartment view and a 4 Mb TAD view want different bin sizes and
# different page widths, and they are the same figure otherwise -- so the knobs
# that legitimately change per region live on the region itself.
loci_raw <- if (length(args) >= 2) {
  setNames(list(args[2]), args[3] %||% gsub("[:-]", "_", args[2]))
} else {
  g$loci
}
if (!length(loci_raw)) {
  stop("no region given: pass one on the command line or fill global.loci in ", cfg_path,
       call. = FALSE)
}

# Keys a locus block interprets itself. Everything else in the block is taken
# as a PANEL override applied to every heatmap of that figure -- which is how a
# region gets its own colour range, matrix or panel heights without editing the
# panels. Both lists are checked, so a typo lands in neither and is reported
# rather than silently ignored.
# Matrix types. Those in COMPARISON_MATRICES need a second file to compare to.
COMPARISON_MATRICES <- c("subtraction", "ratio")
VALID_MATRICES      <- c("counts", "log2", COMPARISON_MATRICES)

LOCUS_KEYS <- c("enabled", "region", "locus", "resolution", "padding_fraction",
                "plot_width", "plot_height", "formato", "dpi", "dir")
# Deliberately NOT overridable per locus: `hic`, `vs`, `name`, `enabled` and
# `preset` identify WHICH panel a block is, so setting them figure-wide would
# point every heatmap at the same file, or name them all alike.
LOCUS_PANEL_KEYS <- c("matrix", "normalize", "title", "height", "min", "max",
                      "colors", "midpoint", "oob", "na_color", "legend")

# A block of panel settings may be written flat, or split into sub-blocks named
# after a matrix type:
#
#   { min: q0.01, max: 12,                     # every panel
#     subtraction: { min: q0.01, max: q0.999 },  # only subtraction panels
#     ratio:       { min: 0.5,   max: 2 } }      # only ratio panels
#
# Contact counts, a log2 ratio, a rescaled difference and a plain ratio live on
# wildly different numeric ranges, so one min/max for all of them is rarely
# what you want -- and repeating the split on every panel is what this avoids.
split_by_matrix <- function(block, what) {
  by_mat <- block[intersect(names(block), VALID_MATRICES)]
  flat   <- block[setdiff(names(block), VALID_MATRICES)]
  # `matrix:` inside a per-matrix block would be circular, so it is not allowed.
  allowed <- setdiff(LOCUS_PANEL_KEYS, "matrix")
  for (mt in names(by_mat)) {
    b <- by_mat[[mt]]
    if (!is.list(b)) {
      stop(what, ".", mt, ": expected a block such as { min: q0.01, max: 12 }", call. = FALSE)
    }
    bad <- setdiff(names(b), allowed)
    if (length(bad)) {
      stop(what, ".", mt, ": unknown key(s) '", paste(bad, collapse = "', '"),
           "'. Allowed here: ", paste(allowed, collapse = ", "), call. = FALSE)
    }
  }
  list(flat = flat, by_matrix = by_mat)
}

loci <- lapply(names(loci_raw), function(nm) {
  e    <- loci_raw[[nm]]
  spec <- if (is.list(e)) e else list(region = e)
  reg  <- spec$region %||% spec$locus
  if (is.null(reg)) {
    stop("global.loci.", nm, ": no region. Write either '", nm, ": chr2:1-2' or '",
         nm, ": { region: chr2:1-2, resolution: 10000 }'", call. = FALSE)
  }
  unknown <- setdiff(names(spec), c(LOCUS_KEYS, LOCUS_PANEL_KEYS, VALID_MATRICES))
  if (length(unknown)) {
    stop("global.loci.", nm, ": unknown key(s) '", paste(unknown, collapse = "', '"),
         "'.\n  Locus keys: ", paste(LOCUS_KEYS, collapse = ", "),
         "\n  Panel overrides allowed here: ", paste(LOCUS_PANEL_KEYS, collapse = ", "),
         "\n  Or a per-matrix block: ", paste(VALID_MATRICES, collapse = ", "),
         call. = FALSE)
  }
  panel_over <- split_by_matrix(spec[setdiff(names(spec), LOCUS_KEYS)],
                                paste0("global.loci.", nm))
  fmt <- tolower(as.character(spec$formato %||% plot_format))
  if (!fmt %in% c("png", "pdf")) {
    stop("global.loci.", nm, ": formato must be 'png' or 'pdf', got '", fmt, "'", call. = FALSE)
  }
  list(
    name    = nm,
    region  = as.character(reg),
    # NA here means "no per-locus override"; the .hic files then fall back to
    # global.resolution, and each file's own `resolution` still wins over both.
    resolution       = num_or_na(spec$resolution),
    padding_fraction = as.numeric(spec$padding_fraction %||% padding_fraction),
    dir     = as.character(spec$dir %||% dir_plots),
    formato = fmt,
    width   = as.numeric(spec$plot_width  %||% plot_width),
    height  = as.numeric(spec$plot_height %||% plot_height),
    dpi     = as.numeric(spec$dpi %||% plot_res),
    # `enabled: false` parks a region without deleting it, so a rerun can
    # regenerate just the figures that changed
    enabled = as_flag(spec$enabled, paste0("global.loci.", nm)),
    # panel settings written on this locus, applied to every heatmap here
    panel_overrides = panel_over)
})

# A CLI region is always drawn: asking for it explicitly outranks the file.
if (length(args) < 2) {
  on_ <- vapply(loci, function(l) l$enabled, logical(1))
  if (!any(on_)) {
    stop("every entry of global.loci has 'enabled: false' in ", cfg_path,
         ": at least one region must be on", call. = FALSE)
  }
  skipped <- vapply(loci[!on_], function(l) l$name, character(1))
  if (length(skipped)) {
    message("skipping disabled loci: ", paste(skipped, collapse = ", "))
  }
  loci <- loci[on_]
}

# ---- hic: named registry of files ------------------------------------------
# Each entry is EITHER a bare path (inheriting the section defaults) OR a map
# that overrides any setting for that file alone:
#
#   files:
#     ysEC: /path/a.hic                                    # inherits everything
#     mESC: { file: /path/b.hic, resolution: 5000, norm: SCALE }
#
# Per-file settings matter because .hic files rarely agree: a deeply sequenced
# library may carry 1 kb bins while a shallower one stops at 10 kb, and not
# every file is balanced with the same algorithm.
hic_raw <- cfg$hic$files
if (!length(hic_raw)) stop("hic.files is empty in ", cfg_path, call. = FALSE)
hic_defaults <- cfg$hic$defaults %||% list()

hic_registry <- setNames(lapply(names(hic_raw), function(nm) {
  e    <- hic_raw[[nm]]
  spec <- if (is.list(e)) e else list(file = e)
  # NA = this file pinned no resolution of its own, so it follows whatever the
  # locus (or the global default) asks for. Precedence, narrowest first:
  #   file's own resolution  >  locus resolution  >  global.resolution
  list(
    name       = nm,
    file       = need_file(spec$file, paste0("hic.files.", nm)),
    resolution = num_or_na(spec$resolution %||% hic_defaults$resolution),
    norm       = as.character(spec$norm %||% hic_defaults$norm %||% cfg$hic$norm %||% "KR"))
}), names(hic_raw))

# Bin size this file will actually be read at, for a given locus.
hic_resolution <- function(hic_name, locus_res) {
  r <- hic_registry[[hic_name]]$resolution
  if (is.na(r)) locus_res else r
}
hic_files <- lapply(hic_registry, `[[`, "file")   # name -> path, for validation

# ---- hic panels: declaration order = top to bottom -------------------------
legend_defaults <- cfg$hic$legend_defaults %||% list()
# Figure-wide panel settings: anything valid on a panel can go here and applies
# to EVERY panel, with each panel still free to override it. This is the place
# for the settings that are usually the same across a whole figure -- matrix,
# normalize, min, max -- so they are written once instead of per panel.
# It also accepts per-matrix sub-blocks, so the four matrix types can carry the
# ranges they actually need without repeating them on every panel:
#   panel_defaults: { normalize: true,
#                     log2:        { min: q0.01, max: 12 },
#                     subtraction: { min: q0.01, max: q0.999 } }
# `hic.defaults` is accepted as the older name.
panel_defaults <- split_by_matrix(cfg$hic$panel_defaults %||% cfg$hic$defaults %||% list(),
                                  "hic.panel_defaults")
panels_cfg <- cfg$hic$panels
if (!length(panels_cfg)) stop("hic.panels is empty in ", cfg_path, call. = FALSE)

# `enabled: false` switches a panel off without deleting it, so one config can
# hold every panel you ever use and each figure just turns on the ones it
# needs. A disabled panel is skipped entirely -- it is never validated and its
# file is never opened -- and the remaining panels close ranks.
panel_presets <- cfg$hic$presets %||% list()
panel_on <- vapply(seq_along(panels_cfg), function(i) {
  what <- paste0("hic.panels[[", i, "]]")
  ef <- merge_settings(panels_cfg[[i]], panel_defaults$flat, panel_presets, what)
  as_flag(ef$enabled, what)
}, logical(1))
panels_cfg <- panels_cfg[panel_on]
if (!length(panels_cfg)) {
  stop("every entry of hic.panels has 'enabled: false' in ", cfg_path,
       ": at least one panel must be on", call. = FALSE)
}

# Default legend title per matrix type; a panel may override it with `title:`.
default_panel_title <- function(name, mtype) {
  switch(mtype,
    counts      = paste0(name, "\nInteractions"),
    log2        = paste0(name, "\nInteractions \n(log2)"),
    subtraction = paste0("subtraction\n", name),
    ratio       = paste0("Ratio\n", name))
}

# Default gradient per matrix type: sequential for a single map, diverging for
# a comparison (which has a meaningful neutral point: 0 for a difference,
# 1 for a ratio).
default_panel_colors <- function(mtype) {
  if (mtype %in% COMPARISON_MATRICES) c("blue", "white", "red")
  else c("white", "orange", "darkred", "black")
}

# Build the panel list for one figure.
#
# `locus_over` are the panel settings written on the locus itself, which is how
# a region gets its own colour range without touching the panels: a 30 Mb
# compartment view and a 4 Mb TAD view rarely share a scale.
#
# Settings resolve in six layers, each narrower than the last:
#
#   hic.panel_defaults  (flat)
#   hic.panel_defaults[matrix]        per matrix type, figure-wide
#   hic.presets[preset]
#   the locus           (flat)
#   the locus[matrix]                 per matrix type, for this region
#   the panel                         always the last word
#
# The per-matrix layers exist because contact counts, a log2 ratio, a rescaled
# difference and a plain ratio are on completely different numeric ranges, so
# a single min/max for a whole figure or region is rarely meaningful.
#
# The panel stays last on purpose. A locus override is about the REGION, so it
# must not silently retype a panel: if the locus said `matrix: log2`, a
# subtraction panel that pinned its own matrix keeps it.
build_panels <- function(locus_over = list(flat = list(), by_matrix = list())) {
  lapply(seq_along(panels_cfg), function(i) {
  what <- paste0("hic.panels[[", i, "]]")
  # Flat layers first, which is enough to know WHICH matrix this panel draws...
  base  <- modifyList(panel_defaults$flat, locus_over$flat)
  mtype <- as.character(
    merge_settings(panels_cfg[[i]], base, panel_presets, what)$matrix %||% "log2")
  # ...then fold in the blocks written for that matrix type, figure-wide first
  # and the locus second, and re-merge so the panel still has the last word.
  base <- modifyList(base, panel_defaults$by_matrix[[mtype]] %||% list())
  base <- modifyList(base, locus_over$by_matrix[[mtype]] %||% list())
  p     <- merge_settings(panels_cfg[[i]], base, panel_presets, what)
  nm    <- as.character(p$name %||% paste0("panel", i))

  if (!mtype %in% VALID_MATRICES) {
    stop("hic.panels[[", i, "]] (", nm, "): unknown matrix '", mtype, "'. Valid: ",
         paste(VALID_MATRICES, collapse = ", "), call. = FALSE)
  }
  hic_a <- as.character(p$hic %||% "")
  if (!nzchar(hic_a) || is.null(hic_files[[hic_a]])) {
    stop("hic.panels[[", i, "]] (", nm, "): 'hic: ", hic_a, "' is not a key of hic.files. ",
         "Known: ", paste(names(hic_files), collapse = ", "), call. = FALSE)
  }
  hic_b <- if (!is.null(p$vs)) as.character(p$vs) else NA_character_
  if (mtype %in% COMPARISON_MATRICES && (is.na(hic_b) || is.null(hic_files[[hic_b]]))) {
    stop("hic.panels[[", i, "]] (", nm, "): matrix '", mtype, "' needs 'vs: <name>' ",
         "pointing at another hic.files key. Known: ",
         paste(names(hic_files), collapse = ", "), call. = FALSE)
  }

  lg  <- p$legend %||% list()
  pos <- unlist(opt(lg, legend_defaults, "position", c(-0.08, 0.5)))

  list(
    name   = nm,
    index  = i,
    hic    = hic_a,
    vs     = hic_b,
    matrix = mtype,
    # Panels sharing a `normalize` group get their colour limits from the
    # POOLED values of the whole group, which is what makes their heatmaps
    # directly comparable. `true` joins a default group; a string names one,
    # so several independent groups can coexist in one figure.
    normalize = if (is.null(p$normalize)) NA_character_
                else if (isTRUE(p$normalize)) "default"
                else if (isFALSE(p$normalize)) NA_character_
                else as.character(p$normalize),
    title  = as.character(p$title %||% default_panel_title(nm, mtype)),
    # relative height of this panel in the final figure
    height = as.numeric(opt(p, list(), "height", 1)),
    # NA = data range, a number = that value, qN = the Nth quantile
    min_spec = parse_limit(p$min, paste0("hic.panels[[", i, "]] (", nm, ") min")),
    max_spec = parse_limit(p$max, paste0("hic.panels[[", i, "]] (", nm, ") max")),
    colors   = as.character(unlist(opt(p, list(), "colors", default_panel_colors(mtype)))),
    # neutral point of a diverging scale: 0 for a difference, 1 for a ratio
    midpoint = as.numeric(opt(p, list(), "midpoint", if (mtype == "ratio") 1 else 0)),
    # What happens to a bin whose value falls outside [min, max]:
    #   censor (default)  it is dropped and drawn in `na_color`. The blank or
    #                     grey cells make the limits of the scale visible, so
    #                     you can see how much of the data you clipped.
    #   squish            it is collapsed onto the nearest end of the scale:
    #                     anything above max takes the max colour, anything
    #                     below min the min colour. Nothing is left uncoloured.
    # `clamp` is accepted as a synonym of squish, `grey`/`na` of censor.
    oob = local({
      v <- tolower(as.character(opt(p, list(), "oob", "censor")))
      switch(v,
        squish = , clamp = "squish",
        censor = , grey = , gray = , na = "censor",
        stop("hic.panels[[", i, "]] (", nm, "): oob must be 'censor' or 'squish', got '",
             v, "'", call. = FALSE))
    }),
    na_color = as.character(opt(p, list(), "na_color",
                                if (mtype %in% COMPARISON_MATRICES) "grey50" else "transparent")),
    show_legend       = as.logical(opt(lg, legend_defaults, "show", TRUE)),
    legend_position   = c(as.numeric(pos[1]), as.numeric(pos[2])),
    legend_direction  = as.character(opt(lg, legend_defaults, "direction", "vertical")),
    legend_title_size = as.numeric(opt(lg, legend_defaults, "title_size", 5)),
    legend_text_size  = as.numeric(opt(lg, legend_defaults, "text_size", 5)),
    legend_key_height = as.numeric(opt(lg, legend_defaults, "key_height", 7)),
    legend_key_width  = as.numeric(opt(lg, legend_defaults, "key_width", 5)),
    # How much room the colourbar is allowed to take. `n_breaks` trims the tick
    # labels, `margin` the white space around the whole legend, and
    # `title_lineheight` the gap between the lines of a multi-line title.
    legend_n_breaks   = as.numeric(opt(lg, legend_defaults, "n_breaks", 0)),
    legend_margin     = as.numeric(unlist(opt(lg, legend_defaults, "margin", 0))),
    legend_title_lineheight = as.numeric(opt(lg, legend_defaults, "title_lineheight", 1)),
    legend_ticks      = as.logical(opt(lg, legend_defaults, "ticks", TRUE)))
})
}

# =============================================================================
# TRACKS: ONE ORDERED LIST, AUTOMATIC VERTICAL LAYOUT
#
# Each entry declares its `type` and its `height`. The band each track occupies
# is computed from the declaration order, the heights and `spacing`, so there
# are no y coordinates in the config to keep consistent.
#
#   band of track i = [baseline_i, baseline_i + height_i]
#   first entry sits at the TOP, each subsequent one `spacing` below the last.
# =============================================================================
tracks_cfg_root <- cfg$tracks %||% list()
track_defaults  <- tracks_cfg_root$defaults %||% list()
track_presets   <- tracks_cfg_root$presets  %||% list()
track_spacing   <- as.numeric(tracks_cfg_root$spacing %||% 0.25)
tracks_label_size <- as.numeric(tracks_cfg_root$label_size %||% 6)
tracks_panel_height <- as.numeric(tracks_cfg_root$panel_height %||% 1)

track_list <- tracks_cfg_root$list
if (!length(track_list)) {
  stop("tracks.list is empty in ", cfg_path,
       ": add at least one entry with 'type: mark | structure | annotation'", call. = FALSE)
}

# `enabled: false` switches a track off without deleting it. Keep every track
# you ever use in the file and flip a couple of booleans per figure. A disabled
# track is skipped before anything is validated, so it may even point at a file
# that is not on this machine, and the tracks below it move up to close the gap.
# `enabled` is resolved through the same three layers, so it can be switched off
# for one track, for a whole preset, or for an entire type at once.
track_on <- vapply(seq_along(track_list), function(i) {
  e    <- track_list[[i]]
  what <- paste0("tracks.list.", names(track_list)[i])
  d    <- track_defaults[[as.character(e$type %||% "")]] %||% list()
  ef   <- merge_settings(e, d, track_presets, what)
  as_flag(ef$enabled, what)
}, logical(1))
track_list <- track_list[track_on]
if (!length(track_list)) {
  stop("every entry of tracks.list has 'enabled: false' in ", cfg_path,
       ": at least one track must be on", call. = FALSE)
}

VALID_TRACK_TYPES <- c("mark", "structure", "annotation")

# --- parse each entry against the defaults of its own type -------------------
tracks <- lapply(names(track_list), function(tn) {
  raw <- track_list[[tn]]
  ty  <- as.character(raw$type %||% "")
  if (!ty %in% VALID_TRACK_TYPES) {
    stop("tracks.list.", tn, ": 'type: ", ty, "' is not valid. Use one of: ",
         paste(VALID_TRACK_TYPES, collapse = ", "), call. = FALSE)
  }
  # Flatten type defaults < preset < entry into one settings list, so every
  # lookup below is a plain read with a hardcoded fallback.
  e <- merge_settings(raw, track_defaults[[ty]] %||% list(), track_presets,
                      paste0("tracks.list.", tn))
  d <- list()

  base <- list(
    name   = tn,
    type   = ty,
    file   = need_file(opt(e, d, "file"), paste0("tracks.list.", tn)),
    height = as.numeric(opt(e, d, "height", 1)),
    # extra gap above this track, on top of the global spacing
    space_above = as.numeric(opt(e, d, "space_above", 0)),
    # label shown on the y axis; defaults to the entry's own name
    label  = as.character(opt(e, d, "label", tn)))

  if (ty == "mark") {
    # signal_scale: "auto" fills the band exactly; a number multiplies
    # log10(1 + signal) directly, in band-height units.
    ss <- opt(e, d, "signal_scale", "auto")
    ss <- if (is.character(ss) && tolower(ss) == "auto") NA_real_ else as.numeric(ss)
    c(base, list(
      color      = as.character(opt(e, d, "color", "grey50")),
      alpha      = as.numeric(opt(e, d, "alpha", 0.9)),
      line_color = as.character(opt(e, d, "line_color", "black")),
      line_width = as.numeric(opt(e, d, "line_width", 0.2)),
      signal_scale = ss,
      n_bins     = as.numeric(opt(e, d, "n_bins", 400)),
      crop       = as.logical(opt(e, d, "crop", TRUE)),
      show_scale = as.logical(opt(e, d, "show_scale", TRUE))))

  } else if (ty == "structure") {
    c(base, list(
      # arc    = semicircle joining the two anchor midpoints (loops)
      # domain = filled bar spanning the whole domain        (TADs)
      style   = as.character(opt(e, d, "style", "arc")),
      # bedpe column handling: read the `select` columns and name them `columns`
      columns = as.character(unlist(opt(e, d, "columns",
                                        c("chr1","start1","end1","chr2","start2","end2")))),
      select  = as.numeric(unlist(opt(e, d, "select", 1:6))),
      # TRUE when the bedpe stores "2" and the genome needs "chr2"
      add_chr = as.logical(opt(e, d, "add_chr", TRUE)),
      # fraction of the band a domain bar fills (arcs always span the band)
      feature_height = as.numeric(opt(e, d, "feature_height", 0.6)),
      # HOW NESTED DOMAINS SHARE THE BAND (style: domain only)
      #   single  every bar on one line, drawn on top of each other
      #   packed  rows, no two OVERLAPPING bars on one row (IGV "squished")
      #   nested  row = nesting depth, so a child always sits under its parent
      layout    = as.character(opt(e, d, "layout", "single")),
      max_rows  = as.numeric(opt(e, d, "max_rows", 0)),      # 0 = as many as needed
      # min gap between two bars sharing a row, as a fraction of the window
      row_gap_frac = as.numeric(opt(e, d, "row_gap_frac", 0.005)),
      fill       = as.character(opt(e, d, "fill", "none")),
      alpha      = as.numeric(opt(e, d, "alpha", 0.25)),
      # per-feature alpha: none (every feature takes `alpha`) | depth | score |
      # span. `alpha_range` is [value_lo, value_hi]; write it backwards to
      # invert the direction.
      alpha_by    = as.character(opt(e, d, "alpha_by", "none")),
      alpha_range = as.numeric(unlist(opt(e, d, "alpha_range", c(0.15, 0.85)))),
      line_color = as.character(opt(e, d, "line_color", "black")),
      line_width = as.numeric(opt(e, d, "line_width", 0.4)),
      line_type  = as.character(opt(e, d, "line_type", "solid")),
      arc_resolution = as.numeric(opt(e, d, "arc_resolution", 60)),
      # keep only structures whose anchor span is in [min_span, max_span]; 0 = off
      min_span  = as.numeric(opt(e, d, "min_span", 0)),
      max_span  = as.numeric(opt(e, d, "max_span", 0)),
      # cap how many are drawn per locus; 0 = draw all
      max_n     = as.numeric(opt(e, d, "max_n", 0)),
      rank_by   = as.character(opt(e, d, "rank_by", "length")),
      rank_desc = as.logical(opt(e, d, "rank_desc", FALSE)),
      score_col = opt(e, d, "score_col"),
      seed      = opt(e, d, "seed")))

  } else {
    c(base, list(
      # model = intron line + exon boxes + thinner UTRs (needs a BED12)
      # box   = one rectangle per feature (works with any BED)
      style = as.character(opt(e, d, "style", "box")),
      # half-height of a feature, as a fraction of the band's half-height:
      # 1 fills the band, 0.2 is a thin bar
      feature_height = as.numeric(opt(e, d, "feature_height", 0.25)),
      # NA = take the color from the BED's itemRgb column
      color      = as.character(opt(e, d, "color", NA_character_)),
      line_color = as.character(opt(e, d, "line_color", "black")),
      line_width = as.numeric(opt(e, d, "line_width", 0.2)),
      utr_thinner = as.logical(opt(e, d, "utr_thinner", TRUE)),
      utr_frac    = as.numeric(opt(e, d, "utr_frac", 0.45)),
      # grow features narrower than this fraction of the window so single-bp
      # exons stay visible; position and overall span are preserved
      min_exon_frac = as.numeric(opt(e, d, "min_exon_frac", 0)),
      # chevrons = repeated ">" marks along the body; triangle = one arrowhead
      direction  = as.character(opt(e, d, "direction", "triangle")),
      chevron_n  = as.numeric(opt(e, d, "chevron_n", 6)),
      arrow_frac = as.numeric(opt(e, d, "arrow_frac", 0.006)),
      # How overlapping features are spread vertically inside the band:
      #   single  every feature on one line (overlaps just draw on top)
      #   packed  real rows, no two overlapping features share a row
      #   jitter  random scatter; shows density, but overlaps are luck
      # Writing `jitter:` without a `layout:` implies layout: jitter, so a
      # track that only asks for scatter still gets it.
      layout    = as.character(opt(e, d, "layout",
                                   if (as.numeric(opt(e, d, "jitter", 0)) > 0) "jitter"
                                   else "single")),
      max_rows  = as.numeric(opt(e, d, "max_rows", 0)),   # 0 = as many as needed
      # minimum horizontal gap between two features sharing a row, as a
      # fraction of the visible window; keeps arrowheads and labels apart
      row_gap_frac = as.numeric(opt(e, d, "row_gap_frac", 0.01)),
      # vertical scatter as a fraction of the band's half-height (layout: jitter)
      jitter = as.numeric(opt(e, d, "jitter", 0)),
      seed   = opt(e, d, "seed"),
      labels = opt_block(e, d, "labels")))
  }
})
names(tracks) <- vapply(tracks, function(t) t$name, character(1))

# --- compute each band from the declaration order ----------------------------
# The first entry is the TOP of the panel; y decreases down the list.
#
# Three knobs control the vertical rhythm, and nothing else does:
#
#   height       (per track)  how tall this track's band is. A mark ribbon
#                             fills its band, so this is the visual weight of
#                             the track.
#   spacing      (global)     gap inserted between EVERY consecutive pair.
#                             Lower it to pack all tracks tighter.
#   space_above  (per track,  extra gap on top of `spacing`, for this track
#                 or per type) only. This is how you separate one group from
#                             another -- e.g. push the marks away from the
#                             annotations -- without tightening everything.
#                             On the first track it becomes top padding.
#
# Because a mark ribbon fills its whole band, the white space you see between
# two marks is exactly `spacing` (+ that track's `space_above`).
track_heights <- vapply(tracks, function(t) t$height, numeric(1))
track_gaps    <- vapply(tracks, function(t) t$space_above, numeric(1))
n_tracks      <- length(tracks)
# total height = all bands + one spacing between each consecutive pair + extras
tracks_total  <- sum(track_heights) + track_spacing * max(0, n_tracks - 1) + sum(track_gaps)

cursor <- tracks_total
for (i in seq_len(n_tracks)) {
  if (i > 1) cursor <- cursor - track_spacing
  cursor <- cursor - track_gaps[i]
  tracks[[i]]$ymin   <- cursor - track_heights[i]   # baseline of the band
  tracks[[i]]$ymax   <- cursor                      # top of the band
  tracks[[i]]$ycen   <- cursor - track_heights[i] / 2
  cursor <- cursor - track_heights[i]
}

# y range of the whole track panel: automatic, with an optional manual override
tl <- as.numeric(unlist(cfg$plot$tracks_limits %||% list()))
tracks_ylim <- if (length(tl) == 2 && all(!is.na(tl))) tl else c(0, tracks_total)

# y axis breaks: one label per track, at the centre of its band
track_axis_dt <- data.table(
  y     = vapply(tracks, function(t) t$ycen, numeric(1)),
  label = vapply(tracks, function(t) t$label, character(1)))

shared_groups <- cfg$shared_yaxis_groups

# The genomic axis drawn on top of the first heatmap: `divide` sets the unit
# (1e6 = Mb), `text_size` the coordinate labels and `title_size` the chromosome
# name above them.
x_axis_cfg       <- cfg$plot$x_axis %||% list()
x_axis_divide    <- as.numeric(x_axis_cfg$divide %||% 1e6)
x_axis_text_size <- as.numeric(x_axis_cfg$text_size  %||% 6)
x_axis_title_size <- as.numeric(x_axis_cfg$title_size %||% 6)

# How much of the rotated square to keep. crop.x trims the left/right corners,
# whose bins pair with positions outside the window; crop.y sets the height of
# the kept slice. Same meaning as the quantiles the original code used.
crop_x <- as.numeric(unlist(cfg$plot$crop$x %||% c(0.2, 0.8)))
crop_y <- as.numeric(cfg$plot$crop$y %||% 0.7)

# Relative heights of the stacked panels: one per heatmap, then the track
# panel. Computed per figure, since a locus may override a panel's height.
panel_proportions <- function(panels) {
  c(vapply(panels, function(p) p$height, numeric(1)), tracks_panel_height)
}

# Convenience views over the track list, by type.
tracks_of <- function(ty) Filter(function(t) t$type == ty, tracks)

# =============================================================================
# HI-C: LOADING AND PANEL MATRICES
# =============================================================================

# A .hic file may name chromosomes "2" or "chr2"; use whichever this file has.
hic_chr_name <- function(path, chr_id) {
  chroms <- strawr::readHicChroms(path)$name
  cand <- c(chr_id, paste0("chr", chr_id), sub("^chr", "", chr_id))
  hit <- cand[cand %in% chroms]
  if (!length(hit)) {
    stop("chromosome '", chr_id, "' not found in ", basename(path),
         ". Available: ", paste(utils::head(chroms, 30), collapse = ", "), call. = FALSE)
  }
  hit[1]
}

# Fail early and clearly if the file cannot serve what its registry entry asks.
check_hic <- function(path, norm, res_bp) {
  norms <- strawr::readHicNormTypes(path)
  if (!norm %in% norms) {
    stop(basename(path), " has no '", norm, "' normalization. Available: ",
         paste(norms, collapse = ", "), " -- set 'norm' for this file in hic.files",
         call. = FALSE)
  }
  res <- strawr::readHicBpResolutions(path)
  if (!res_bp %in% res) {
    stop(basename(path), " has no ", res_bp, " bp resolution. Available: ",
         paste(res, collapse = ", "), " -- set 'resolution' for this file in hic.files",
         call. = FALSE)
  }
  invisible(TRUE)
}

# Whole-chromosome matrix for ONE hic file, at that file's own resolution and
# normalisation. Cached per (file, chromosome, resolution, norm), so several
# loci on the same chromosome -- and several panels sharing a file -- pay the
# straw cost once.
hic_cache <- new.env(parent = emptyenv())
load_hic_chr <- function(hic_name, chr_id, res_bp) {
  h   <- hic_registry[[hic_name]]
  key <- paste(hic_name, chr_id, res_bp, h$norm, sep = "|")
  if (!is.null(hic_cache[[key]])) return(hic_cache[[key]])

  path <- h$file
  size <- chr.list[chr == chr_id, size]
  if (!length(size)) {
    stop("chromosome '", chr_id, "' is not in genome '", genome_name, "'", call. = FALSE)
  }
  # Expected number of populated bins in the upper triangle of the chromosome.
  # Computed at the resolution actually used, so the mean below is comparable
  # across files even when their bin sizes differ.
  bins.chr <- (round(size / res_bp)^2) / 2

  check_hic(path, h$norm, res_bp)
  nm <- hic_chr_name(path, chr_id)
  message("  loading ", hic_name, " (chr", chr_id, ", ", res_bp, " bp, ", h$norm, ")")
  dt <- strawr::straw(h$norm, path, nm, nm, "BP", res_bp, matrix = "observed")
  setDT(dt)

  avg <- sum(dt$counts, na.rm = TRUE) / bins.chr    # mean contacts per bin pair
  dt[, counts := counts / avg]                      # depth-independent units
  dt[, log2 := log2(counts)]
  setkey(dt, x, y)

  res <- list(dt = dt, avg = avg)
  hic_cache[[key]] <- res
  res
}

# Values for one panel over the locus grid, as x / y / value.
# Bins with no contacts are filled with the matrix type's neutral value, so an
# empty bin reads as "nothing special here" rather than as missing data.
panel_locus_map <- function(panel, chr_id, win_start, win_end, res_bp) {
  # Each panel tiles the shared window at its OWN resolution.
  anchors <- seq.int(win_start, win_end, res_bp)
  grid <- CJ(x = anchors, y = anchors, sorted = FALSE)
  setkey(grid, x, y)

  a <- load_hic_chr(panel$hic, chr_id, res_bp)

  if (panel$matrix %in% c("counts", "log2")) {
    m <- a$dt[grid]
    m[, value := get(panel$matrix)]
    neutral <- 0
  } else {
    b <- load_hic_chr(panel$vs, chr_id, res_bp)
    m <- merge(a$dt, b$dt, by = c("x", "y"), all = TRUE, suffixes = c("", "_vs"))
    setkey(m, x, y)
    m <- m[grid]
    if (panel$matrix == "subtraction") {
      # Rescale the difference of two normalised maps back to a contact-like
      # magnitude, using the mean of both files' densities.
      avg_sum <- (a$avg / 2) + (b$avg / 2)
      m[, value := (counts - counts_vs) * avg_sum]
      neutral <- 0
    } else {
      m[, value := counts / counts_vs]
      neutral <- 1     # a ratio of 1 is "no change"
    }
  }
  m[, value := nafill(value, fill = neutral)]
  m[, .(x, y, value)]
}

# =============================================================================
# HEATMAP RENDERING
# =============================================================================

# The four corners of one matrix bin, rotated -45 degrees so the diagonal of
# the contact matrix becomes the horizontal genomic axis.
to_polygon <- function(x, y, counts, resol) {
  half <- resol / 2
  cx <- c(x - half, x + half, x + half, x - half)
  cy <- c(y - half, y - half, y + half, y + half)
  data.frame(u = (cx + cy) / sqrt(2), v = (cy - cx) / sqrt(2),
             counts = counts, id = paste(x, y, sep = "_"), corner = 1:4)
}

# Build the fill scale for a panel. A comparison uses a diverging ramp pinned
# at its neutral point (0 or 1); a single map uses a sequential ramp.
# `scale_vals` is what the limits are computed FROM. For a normal panel it is
# that panel's own values; for a panel in a `normalize` group it is the pooled
# values of the whole group, so every member -- given the same min/max spec --
# lands on identical limits and the heatmaps become directly comparable.
# Quantiles follow the same rule: q0.99 of the pool, not of one sample.
panel_scale <- function(panel, scale_vals) {
  lo <- resolve_limit(panel$min_spec, scale_vals, min)
  hi <- resolve_limit(panel$max_spec, scale_vals, max)
  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
    # Degenerate range (an empty or constant window); fall back to the data.
    lo <- min(scale_vals, na.rm = TRUE); hi <- max(scale_vals, na.rm = TRUE)
    if (!is.finite(lo) || !is.finite(hi)) { lo <- 0; hi <- 1 }
    if (hi <= lo) hi <- lo + 1
  }
  # squish collapses out-of-range values onto the nearest limit; censor drops
  # them so they render in na_color, which keeps the clipping visible.
  oob_fun <- if (identical(panel$oob, "squish")) scales::squish else scales::censor

  # Fewer ticks make the colourbar shorter, which is the main lever when
  # neighbouring panels' legends start colliding.
  brks <- if (panel$legend_n_breaks >= 2) {
    scales::breaks_pretty(n = panel$legend_n_breaks)(c(lo, hi))
  } else waiver()

  args <- list(name = panel$title, colours = panel$colors, limits = c(lo, hi),
               breaks = brks, oob = oob_fun, na.value = panel$na_color)
  if (panel$matrix %in% COMPARISON_MATRICES) {
    # `values` places the middle color exactly at the neutral point, so
    # asymmetric limits do not shift where "no change" is drawn.
    args$values <- scales::rescale(c(lo, panel$midpoint, hi))
  }
  do.call(scale_fill_gradientn, args)
}

# One rotated heatmap, cropped to limits supplied by the caller.
#
# The crop is computed once per locus from the window itself (see crop_limits),
# NOT from this panel's own data. That is what lets panels of different
# resolutions line up exactly: a quantile taken over each panel's own bins
# would land half a bin apart for a 5 kb panel and a 20 kb one.
matriz_horizontal <- function(panel, plot_df, region, resol, lims, scale_vals) {
  plot_df <- plot_df[, .(x = as.numeric(x), y = as.numeric(y), counts = value)]

  # Only the upper triangle is drawn; the contact matrix is symmetric.
  polygonos <- plot_df[y >= x][, to_polygon(x, y, counts, resol), by = .(x, y)]

  p <- ggplot(polygonos, aes(u, v, group = id, fill = counts)) +
    geom_polygon() +
    coord_cartesian(clip = "off") +
    panel_scale(panel, scale_vals) +
    ggtitle(region) +
    scale_y_continuous(limits = c(0, lims$v_max), expand = c(0, 0)) +
    scale_x_continuous(limits = c(lims$u_min, lims$u_max),
                       position = "top", expand = c(0, 0),
                       # u is a rotated coordinate: /sqrt(2) turns it back into bp
                       labels = function(u) round(u / sqrt(2) / x_axis_divide, 2)) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = x_axis_text_size, angle = 0, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.text.y = element_blank(), axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = x_axis_title_size),
      axis.line = element_blank(),
      axis.ticks.x = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length.x = unit(3, "pt"),
      axis.ticks.y = element_blank(), axis.ticks.length.y = element_blank(),
      panel.border = element_rect(linetype = "solid", fill = NA),
      plot.margin = margin(0, 0, 0, 20),
      legend.position = if (panel$show_legend) panel$legend_position else "none",
      legend.direction = panel$legend_direction,
      legend.title = element_text(size = panel$legend_title_size, face = "bold",
                                  hjust = 0.5, lineheight = panel$legend_title_lineheight),
      legend.text  = element_text(size = panel$legend_text_size,  face = "bold", hjust = 0.5),
      legend.key.height = unit(panel$legend_key_height, "pt"),
      legend.key.width  = unit(panel$legend_key_width, "pt"),
      # recycled to the 4 sides, so `margin: 0` or `margin: [0,2,0,2]` both work
      legend.margin = do.call(margin, as.list(rep_len(panel$legend_margin, 4))),
      legend.box.spacing = unit(0, "pt")) +
    guides(fill = guide_colourbar(ticks = panel$legend_ticks,
                                  draw.ulim = FALSE, draw.llim = FALSE))
  p
}

# Crop limits for a locus, derived analytically from the window rather than
# from any panel's data.
#
# On a uniform grid spanning [ws, we], the rotated coordinates are uniform too:
#   u = (x + y)/sqrt(2)  is uniform over [2*ws, 2*we]/sqrt(2)
#   v = (y - x)/sqrt(2)  is uniform over [-W, W]/sqrt(2),  W = we - ws
# so the q-th quantile of each has a closed form that does not depend on the
# bin size at all. `crop.x` keeps the middle of the map (the corners hold bins
# whose partner falls outside the window); `crop.y` sets how tall a slice of
# the rotated square is kept.
crop_limits <- function(win_start, win_end, crop_x, crop_y) {
  W <- win_end - win_start
  vis_start <- win_start + crop_x[1] * W
  vis_end   <- win_start + crop_x[2] * W
  list(u_min = sqrt(2) * vis_start,
       u_max = sqrt(2) * vis_end,
       v_max = (2 * crop_y - 1) * W / sqrt(2),
       visible = list(start = floor(vis_start), end = ceiling(vis_end)))
}

# =============================================================================
# MARKS: BIGWIG SIGNAL
# =============================================================================

# Mean score of the overlapping bigwig intervals per bin, one row per bin.
bin_marks <- function(region_gr) {
  mk_tracks <- tracks_of("mark")
  if (!length(mk_tracks)) return(data.table())
  rbindlist(lapply(mk_tracks, function(t) {
    bins <- tile(region_gr, n = t$n_bins)[[1]]
    bw   <- import(t$file, which = region_gr)
    sig  <- numeric(length(bins))
    if (length(bw)) {
      ov <- findOverlaps(bins, bw)
      if (length(ov)) {
        agg <- tapply(bw$score[subjectHits(ov)], queryHits(ov), mean)
        sig[as.integer(names(agg))] <- agg
      }
    }
    data.table(start = start(bins), end = end(bins), signal = sig, mark = t$name)
  }))
}

# =============================================================================
# STRUCTURES: TADs AND LOOPS FROM BEDPE
# =============================================================================

# Read a bedpe and reduce it to the pieces the plot needs: `start`/`end` are the
# MIDPOINTS of the two anchors (the points a loop actually connects), while
# `dstart`/`dend` are the outer edges (the extent a domain covers).
load_structure <- function(s) {
  dt <- fread(s$file, select = s$select, col.names = s$columns)
  chr <- if (isTRUE(s$add_chr)) paste0("chr", dt$chr1) else as.character(dt$chr1)
  sc  <- if (!is.null(s$score_col) && s$score_col %in% names(dt)) as.numeric(dt[[s$score_col]])
         else NA_real_
  data.table(chr = chr,
             start  = (dt$start1 + dt$end1) / 2,
             end    = (dt$start2 + dt$end2) / 2,
             dstart = pmin(dt$start1, dt$start2),
             dend   = pmax(dt$end1, dt$end2),
             score  = sc)
}

# NESTING DEPTH: how many other structures of the same set contain this one.
# 0 = outermost, 1 = inside one domain, and so on. This is what tells a
# sub-TAD from the TAD around it, and it is the value `alpha_by: depth` and
# `layout: nested` are built on.
#
# Containment is counted directly rather than walked as a tree, because TAD
# callers happily emit domains that partially overlap: those are not a
# hierarchy at all, and a stack-based descent would have to invent a parent for
# them. Identical coordinates are broken by row order, so two copies of one
# domain never each count as the other's parent.
nesting_depth <- function(lo, hi) {
  n <- length(lo)
  if (!n) return(integer(0))
  span <- hi - lo
  idx  <- seq_len(n)
  vapply(idx, function(i) {
    sum(lo <= lo[i] & hi >= hi[i] &
        (span > span[i] | (span == span[i] & idx < i)))
  }, integer(1))
}

# Per-feature alpha. `none` gives every feature the track's own `alpha`;
# otherwise the chosen variable is stretched onto `alpha_range`, so the whole
# range is used whatever the locus holds.
#
# `depth` is anchored, not stretched: depth 0 always takes the FIRST end of the
# range, so an outermost TAD looks the same in a locus nested three deep as in
# one nested six deep, and only the levels above it get more intense. For
# `score` and `span` there is no such natural zero, so the visible minimum takes
# that end instead. Either way, write alpha_range backwards to invert it.
structure_alpha <- function(s, sc) {
  by <- tolower(s$alpha_by %||% "none")
  n  <- nrow(sc)
  if (by %in% c("none", "", "na")) return(rep(s$alpha, n))

  rg <- s$alpha_range
  if (length(rg) != 2 || any(is.na(rg))) {
    stop("tracks.list.", s$name, ": 'alpha_range' must be two numbers, [lo, hi]",
         call. = FALSE)
  }
  v <- switch(by,
    depth = sc$depth,
    score = sc$score,
    span  = sc$xhi - sc$xlo,
    stop("tracks.list.", s$name, ": 'alpha_by: ", s$alpha_by,
         "' is not valid. Use one of: none, depth, score, span", call. = FALSE))

  if (all(is.na(v))) {
    warning("tracks.list.", s$name, ": 'alpha_by: ", by,
            "' has no values here; falling back to the track's own alpha", call. = FALSE)
    return(rep(s$alpha, n))
  }
  t <- if (by == "depth") {
         v / max(1, max(v, na.rm = TRUE))            # anchored at depth 0
       } else {
         rng <- range(v, na.rm = TRUE)
         if (diff(rng) <= 0) rep(0, n) else (v - rng[1]) / diff(rng)
       }
  t[is.na(t)] <- 0
  pmin(1, pmax(0, rg[1] + t * (rg[2] - rg[1])))
}

# Structures of one track that fall inside the visible window, after the span
# filters and the max_n cap.
select_structures <- function(s, chr_label, win_start, win_end) {
  sc <- load_structure(s)[chr == chr_label]
  if (!nrow(sc)) return(sc)
  sc <- sc[end > win_start & start < win_end]          # any overlap with the view
  if (!nrow(sc)) return(sc)

  if (s$min_span > 0) sc <- sc[(end - start) >= s$min_span]
  if (s$max_span > 0) sc <- sc[(end - start) <= s$max_span]
  if (!nrow(sc)) return(sc)

  # The extent this track actually draws: the whole domain for a bar, the two
  # anchor midpoints for an arc. Nesting is judged on what the reader sees, so
  # both depth and `alpha_by: span` read the same coordinates as the geometry.
  if (identical(s$style, "domain")) sc[, `:=`(xlo = dstart, xhi = dend)]
  else                              sc[, `:=`(xlo = start,  xhi = end)]

  # Depth is measured BEFORE the cap, on everything that reaches the window: a
  # sub-TAD does not stop being a sub-TAD because `max_n` dropped the domain
  # around it, and letting the cap renumber the levels would make the same
  # region change colour with an unrelated setting. Nothing outside the window
  # can contain something inside it, so this set is already complete.
  sc[, depth := nesting_depth(xlo, xhi)]

  # Cap how many are drawn, keeping the most relevant ones first.
  if (s$max_n > 0 && nrow(sc) > s$max_n) {
    ord <- switch(s$rank_by,
      score  = order(sc$score, decreasing = s$rank_desc, na.last = TRUE),
      random = { if (!is.null(s$seed)) set.seed(s$seed); sample(seq_len(nrow(sc))) },
      order(sc$end - sc$start, decreasing = s$rank_desc))   # "length"
    sc <- sc[utils::head(ord, s$max_n)]
  }
  sc
}

# Geometry for one structure track: semicircular arcs, or domain bars.
build_structure_track <- function(s, chr_label, win_start, win_end) {
  sc <- select_structures(s, chr_label, win_start, win_end)
  if (!nrow(sc)) {
    warning("tracks.list.", s$name, ": nothing to draw in the visible region", call. = FALSE)
    return(list(arcs = data.table(), domains = data.table()))
  }

  alp <- structure_alpha(s, sc)

  # ggplot applies a shape's alpha to its FILL, never to its outline, so on its
  # own `alpha_by` would leave every feature wearing the same solid border and
  # an unfilled track would show no shading at all. The outline colour is
  # therefore built per feature here, with the same alpha baked in -- but only
  # while `alpha_by` is on. With it off, `alpha` is the track's single fill
  # alpha, and pushing that onto the outline as well would wash out every arc
  # track that already relies on a solid line over a translucent fill.
  vary <- !tolower(s$alpha_by %||% "none") %in% c("none", "", "na")
  line <- if (identical(s$line_color, "none")) NA_character_
          else if (vary) scales::alpha(s$line_color, alp)
          else s$line_color

  if (s$style == "domain") {
    # A filled bar spanning the whole domain: the natural TAD representation in
    # a linear browser panel.
    #
    # ROWS. `single` keeps every bar on one line, which is honest only while
    # the domains are disjoint -- give it a nested caller and the sub-TADs
    # vanish under their parents. `packed` and `nested` split the band into
    # rows instead, and, as in the annotation track, feature height is relative
    # to the ROW, so adding rows makes the bars thinner rather than letting
    # them overflow the band.
    n <- nrow(sc)
    rows <- switch(s$layout,
      single = rep(1L, n),
      # IGV's "squished": first free row, left to right. Compact, but a row
      # says nothing about hierarchy -- two bars share one only because they
      # happen not to overlap.
      packed = pack_rows(sc$xlo, sc$xhi, s$row_gap_frac * (win_end - win_start),
                         s$max_rows),
      # Row = nesting depth, so every outermost domain is on row 1 and a child
      # is always exactly one row below its parent. Wastes more space than
      # packing, and reads as a hierarchy because it is one.
      nested = { r <- sc$depth + 1L
                 if (s$max_rows > 0) pmin(r, as.integer(s$max_rows)) else r },
      stop("tracks.list.", s$name, ": 'layout: ", s$layout,
           "' is not valid for a structure track. Use one of: single, packed, nested",
           call. = FALSE))

    n_rows <- max(rows)
    row_h  <- s$height / n_rows
    h      <- row_h / 2 * s$feature_height     # half-height of a drawn bar
    # row 1 at the TOP of the band, like a genome browser
    ycen   <- s$ymax - (rows - 0.5) * row_h

    return(list(
      arcs = data.table(),
      domains = data.table(
        xmin = sc$xlo, xmax = sc$xhi,
        ymin = ycen - h, ymax = ycen + h,
        alpha = alp, line = line, depth = sc$depth,
        id = paste0(s$name, "_", seq_len(n)), track = s$name)))
  }

  # Arc: a semicircle from anchor midpoint to anchor midpoint, rising from the
  # band's baseline. Heights are normalised by the widest arc in view, so the
  # tallest one exactly fills the band and none can escape it.
  #
  # Deepest last, so when the arcs are filled the small nested ones are not
  # buried under the wide one that contains them.
  ord <- order(sc$depth)
  sc <- sc[ord]; alp <- alp[ord]
  if (length(line) > 1) line <- line[ord]
  res <- max(3, s$arc_resolution)
  widest <- max(sc$end - sc$start)
  out <- lapply(seq_len(nrow(sc)), function(j) {
    x1 <- sc$start[j]; x2 <- sc$end[j]
    if (x2 <= x1) return(NULL)
    xs  <- seq(x1, x2, length.out = res)
    rad <- (x2 - x1) / 2; xc <- (x1 + x2) / 2
    h <- sqrt(pmax(0, rad^2 - (xs - xc)^2)) / (widest / 2) * s$height
    data.table(x = xs, y = s$ymin + h,
               alpha = alp[j], line = if (length(line) > 1) line[j] else line,
               depth = sc$depth[j],
               id = paste0(s$name, "_", j), track = s$name)
  })
  list(arcs = rbindlist(Filter(Negate(is.null), out)), domains = data.table())
}

# =============================================================================
# ANNOTATIONS: ANY BED, AS A GENE MODEL OR AS BOXES
# =============================================================================

# Enforce a minimum drawn width so a 1 bp exon is still visible in a 4 Mb
# window. Each short feature grows symmetrically about its midpoint, so its
# position and the overall span of the gene are preserved.
# Assign features to rows so that no two features on the same row overlap --
# the "packed" layout every genome browser uses. Greedy and deterministic:
# walk the features left to right and drop each one on the first row whose last
# feature ended at least `pad` bp earlier.
#
# This is what `jitter` only approximates: jitter scatters at random, so it
# separates features by luck and shuffles on every run, while packing
# guarantees no overlap and is reproducible.
pack_rows <- function(start, end, pad = 0, max_rows = 0) {
  n <- length(start)
  row <- integer(n)
  row_end <- numeric(0)
  for (i in order(start)) {
    placed <- FALSE
    for (r in seq_along(row_end)) {
      if (start[i] > row_end[r] + pad) {
        row[i] <- r; row_end[r] <- end[i]; placed <- TRUE; break
      }
    }
    if (!placed) {
      if (max_rows > 0 && length(row_end) >= max_rows) {
        # Out of rows: put it on whichever row frees up soonest. Overlap is
        # accepted here, but only once the cap has actually been reached.
        r <- which.min(row_end)
        row[i] <- r; row_end[r] <- max(row_end[r], end[i])
      } else {
        row_end <- c(row_end, end[i]); row[i] <- length(row_end)
      }
    }
  }
  row
}

widen <- function(xmin, xmax, minw) {
  if (minw <= 0) return(list(xmin = xmin, xmax = xmax))
  short <- (xmax - xmin) < minw
  mid <- (xmin + xmax) / 2
  xmin[short] <- mid[short] - minw / 2
  xmax[short] <- mid[short] + minw / 2
  list(xmin = xmin, xmax = xmax)
}

# Choose which features get a text label. This is what keeps a dense RefSeq
# track readable: collapse isoforms, drop tiny features, and cap the count.
select_labels <- function(dt, lab, window) {
  if (!isTRUE(lab$enabled %||% FALSE) || !nrow(dt)) return(dt[0])
  d <- copy(dt)
  d[, width := end - start]

  # Names that must always be labelled, whatever the filters say.
  keep_always <- as.character(unlist(lab$include %||% character(0)))
  drop_names  <- as.character(unlist(lab$exclude %||% character(0)))
  if (length(drop_names)) d <- d[!gene %in% drop_names]
  if (!nrow(d)) return(d)

  # Collapse isoforms sharing a name down to the longest one, so a gene is
  # named once instead of once per transcript. This runs BEFORE the include
  # list is split off: `include` means "make sure this gene is labelled", not
  # "label every transcript of it", and a gene with several isoforms would
  # otherwise come back once per isoform.
  if (isTRUE(lab$dedup %||% TRUE)) {
    setorder(d, -width)
    d <- d[, .SD[1], by = gene]
  }

  # Named genes bypass the size filter and the cap below, but not dedup.
  forced <- d[gene %in% keep_always]
  pool   <- d[!gene %in% keep_always]
  # Ignore features too small to be worth naming, in bp or relative to the view.
  min_w <- as.numeric(lab$min_width %||% 0)
  if (min_w > 0) pool <- pool[width >= min_w]
  min_wf <- as.numeric(lab$min_width_frac %||% 0)
  if (min_wf > 0) pool <- pool[width >= min_wf * window]

  # Cap the number of labels, keeping the most representative first.
  max_n <- as.numeric(lab$max_n %||% 0)
  if (max_n > 0 && nrow(pool) > max_n) {
    ord <- switch(as.character(lab$rank_by %||% "length"),
      score    = order(pool$score, decreasing = !isFALSE(lab$rank_desc), na.last = TRUE),
      position = order(pool$start),
      random   = { if (!is.null(lab$seed)) set.seed(lab$seed); sample(seq_len(nrow(pool))) },
      order(pool$width, decreasing = !isFALSE(lab$rank_desc)))   # "length"
    pool <- pool[utils::head(ord, max_n)]
  }
  rbind(forced, pool, fill = TRUE)
}

# Build every geometry a single annotation track needs.
build_annotation_track <- function(a, chr_id, region_gr, window) {
  empty <- data.table()
  gr <- import(a$file, which = region_gr)
  gr <- gr[seqnames(gr) == paste0("chr", chr_id)]
  gr <- subsetByOverlaps(gr, region_gr)
  if (!length(gr)) {
    warning("tracks.list.", a$name, ": no features in the visible region", call. = FALSE)
    return(list(boxes = empty, introns = empty, exons = empty, cds = empty,
                chevrons = empty, triangles = empty, labels = empty,
                label_zone = list(centre = a$ycen, slack = a$height / 2)))
  }

  mc <- mcols(gr)
  # Per-feature color: an explicit track color wins, otherwise the BED's own
  # itemRgb column, otherwise black.
  col <- if (!is.na(a$color)) rep(a$color, length(gr))
         else if (!is.null(mc$itemRgb)) vapply(mc$itemRgb, rgb_from_string, character(1))
         else rep("#000000", length(gr))
  nms  <- if (!is.null(mc$name))  as.character(mc$name) else rep(NA_character_, length(gr))
  scr  <- if (!is.null(mc$score)) as.numeric(mc$score)  else rep(NA_real_, length(gr))
  strd <- as.character(strand(gr))

  # ---- vertical placement inside this track's band --------------------------
  # `packed` splits the band into as many rows as the features need; `single`
  # and `jitter` keep one row. Feature height is always relative to the ROW, so
  # adding rows makes each glyph proportionally thinner instead of overflowing.
  n_feat <- length(gr)
  if (a$layout == "packed") {
    rows   <- pack_rows(start(gr), end(gr), a$row_gap_frac * window, a$max_rows)
    n_rows <- max(rows)
  } else {
    rows   <- rep(1L, n_feat)
    n_rows <- 1L
  }

  # `labels.reserve` takes a slice off the TOP of the band and keeps it clear of
  # features, so the labels have somewhere to live. Without it the only way to
  # get labels off the genes is to nudge them out of the band by hand, and any
  # slack you then give repulsion lets it drop them back onto the features.
  # With it, the separation holds by construction: features pack below, labels
  # sit above, and neither can reach the other.
  lab_cfg <- a$labels
  reserve <- max(0, min(0.9, as.numeric(lab_cfg$reserve %||% 0)))
  feat_h  <- a$height * (1 - reserve)          # height left for the features
  feat_top <- a$ymin + feat_h                  # boundary between the two zones

  row_h  <- feat_h / n_rows                    # height of one row
  half   <- row_h / 2                          # half-height of a row
  h      <- half * a$feature_height            # half-height of a drawn feature
  utrh   <- if (isTRUE(a$utr_thinner)) h * a$utr_frac else h
  minw   <- a$min_exon_frac * window

  # row 1 sits at the TOP of the feature zone, like a genome browser
  ybase <- feat_top - (rows - 0.5) * row_h
  if (a$layout == "jitter" && a$jitter > 0) {
    if (!is.null(a$seed)) set.seed(a$seed)
    ybase <- ybase + runif(n_feat, -half * a$jitter, half * a$jitter)
  }

  # Where the labels are allowed to be, as an interval [zone_lo, zone_hi].
  #
  # Repulsion is free to put a label anywhere in this interval, and it uses
  # that freedom to resolve collisions -- which is why, given the whole strip,
  # some labels end up high and others come to rest right on the features.
  # `gap` raises the FLOOR of the interval, so no label can settle closer than
  # that to the genes no matter how crowded the locus is. `ylim` then shrinks
  # the interval symmetrically, and `offset` slides the whole thing.
  band_half <- a$height / 2
  if (reserve > 0) {
    strip_h <- a$height * reserve
    gap     <- max(0, min(0.95, as.numeric(lab_cfg$gap %||% 0)))
    zone_lo <- feat_top + strip_h * gap
    zone_hi <- a$ymax
  } else {
    zone_lo <- a$ycen - band_half
    zone_hi <- a$ycen + band_half
  }
  zone_half   <- (zone_hi - zone_lo) / 2
  zone_centre <- (zone_lo + zone_hi) / 2 + band_half * as.numeric(lab_cfg$offset %||% 0)
  label_zone <- list(centre = zone_centre,
                     slack  = zone_half * as.numeric(lab_cfg$ylim %||% 1))

  feat <- data.table(start = start(gr), end = end(gr), strand = strd,
                     gene = nms, score = scr, color = col, y = ybase, track = a$name)

  boxes <- introns <- exons <- cds <- chevrons <- triangles <- empty

  if (a$style == "model" && !is.null(mc$blocks)) {
    # --- full gene model: intron line, exon boxes, thicker CDS ---------------
    blocks <- mc$blocks
    thick  <- mc$thick                                # CDS extent, may be absent

    introns <- data.table(x = feat$start, xend = feat$end,
                          y = feat$y, yend = feat$y, track = a$name)

    ex <- rbindlist(lapply(seq_along(gr), function(k) {
      b <- blocks[[k]]
      if (!length(b)) return(NULL)
      # BED block coordinates are relative to the feature start
      xs <- feat$start[k] + start(b) - 1
      xe <- feat$start[k] + end(b) - 1
      w  <- widen(xs, xe, minw)
      data.table(xmin = w$xmin, xmax = w$xmax,
                 ymin = feat$y[k] - utrh, ymax = feat$y[k] + utrh,
                 color = feat$color[k], track = a$name, idx = k)
    }))
    exons <- if (is.null(ex)) empty else ex

    # The coding part of each exon is drawn taller, which is what makes UTRs
    # read as the thinner ends of the transcript.
    if (!is.null(thick) && nrow(exons)) {
      cd <- rbindlist(lapply(seq_along(gr), function(k) {
        e <- exons[idx == k]
        if (!nrow(e)) return(NULL)
        c1 <- pmax(e$xmin, start(thick)[k]); c2 <- pmin(e$xmax, end(thick)[k])
        keep <- c2 > c1
        if (!any(keep)) return(NULL)
        w <- widen(c1[keep], c2[keep], minw)
        data.table(xmin = w$xmin, xmax = w$xmax,
                   ymin = feat$y[k] - h, ymax = feat$y[k] + h,
                   color = feat$color[k], track = a$name)
      }))
      cds <- if (is.null(cd)) empty else cd
    }
  } else {
    # --- plain boxes: works for any BED (peaks, motifs, enhancers, ...) ------
    if (a$style == "model") {
      warning("tracks.list.", a$name,
              ": style 'model' needs a BED12 with block columns; drawing boxes instead",
              call. = FALSE)
    }
    w <- widen(feat$start, feat$end, minw)
    boxes <- data.table(xmin = w$xmin, xmax = w$xmax,
                        ymin = feat$y - h, ymax = feat$y + h,
                        color = feat$color, track = a$name)
  }

  # --- strand indicator ------------------------------------------------------
  if (a$direction == "chevrons") {
    # Repeated ">" marks along the body: readable at any zoom, and they do not
    # add width to the feature the way an arrowhead does.
    n  <- max(1, a$chevron_n)
    d  <- window * 0.004
    hh <- h * 0.7
    ch <- rbindlist(lapply(seq_len(nrow(feat)), function(k) {
      if (!feat$strand[k] %in% c("+", "-")) return(NULL)
      xcv <- seq(feat$start[k], feat$end[k], length.out = n + 2)[2:(n + 1)]
      ax  <- if (feat$strand[k] == "+") xcv - d else xcv + d
      data.table(x = c(ax, ax),
                 y = c(rep(feat$y[k] - hh, n), rep(feat$y[k] + hh, n)),
                 xend = c(xcv, xcv), yend = rep(feat$y[k], 2 * n), track = a$name)
    }))
    chevrons <- if (is.null(ch)) empty else ch
  } else if (a$direction == "triangle") {
    aw <- window * a$arrow_frac
    tr <- rbindlist(lapply(seq_len(nrow(feat)), function(k) {
      if (!feat$strand[k] %in% c("+", "-")) return(NULL)
      tx <- if (feat$strand[k] == "+") c(feat$end[k], feat$end[k] + aw, feat$end[k])
            else                       c(feat$start[k], feat$start[k] - aw, feat$start[k])
      data.table(x = tx, y = c(feat$y[k] - h, feat$y[k], feat$y[k] + h),
                 color = feat$color[k], id = paste0(a$name, "_", k), track = a$name)
    }))
    triangles <- if (is.null(tr)) empty else tr
  }

  # --- labels ----------------------------------------------------------------
  lab <- a$labels
  labels <- empty
  cand <- feat[!is.na(gene) & nzchar(gene)]
  if (nrow(cand)) {
    sel <- select_labels(cand, lab, window)
    if (nrow(sel)) {
      # label jitter and offset are fractions of the band's half-height too
      # `offset` and `jitter` are fractions of the BAND half-height, not of the
      # row. A label is pushed away from the features as a group, and with
      # packed rows a row-relative offset would shrink with every extra row --
      # exactly when you need it most.
      jt <- zone_half * as.numeric(lab$jitter %||% 0)
      if (!is.null(lab$seed)) set.seed(lab$seed)
      labels <- data.table(
        x = (sel$start + sel$end) / 2,
        # The anchor stays ON the feature. ggrepel draws its leader line from
        # the label back to this point, so baking the offset into y here would
        # aim every line at empty space above the gene instead of at the gene.
        y = sel$y,
        # Nudge each label from its feature up to the label zone, so repulsion
        # starts where the label belongs instead of on top of the gene.
        nudge = (label_zone$centre - sel$y) +
                (if (jt > 0) runif(nrow(sel), -jt, jt) else 0),
        label = sel$gene, track = a$name)
    }
  }

  list(boxes = boxes, introns = introns, exons = exons, cds = cds,
       chevrons = chevrons, triangles = triangles, labels = labels,
       # computed once here so the renderer never re-derives the geometry
       label_zone = label_zone)
}

# =============================================================================
# RENDER ONE LOCUS
# =============================================================================
render_locus <- function(locus) {
  nomen <- locus$name
  reg <- parse_region(locus$region)
  chr_id    <- sub("^chr", "", reg$chr)
  chr_label <- paste0("chr", chr_id)
  message("== ", nomen, " (", chr_label, ":", reg$start, "-", reg$end, ") ==")

  # Panels are rebuilt per figure, so the locus can carry panel overrides
  # (min, max, matrix, normalize, colors, height, ...) that apply to all of them.
  panels <- build_panels(locus$panel_overrides)

  # Bin size for this locus: the locus override if it has one, else the global
  # default. A .hic that pinned its own resolution ignores both.
  locus_res <- if (is.na(locus$resolution)) resol else locus$resolution
  panel_res <- vapply(panels, function(p) hic_resolution(p$hic, locus_res), numeric(1))

  # A comparison merges two matrices bin by bin, so both sides must agree.
  for (i in seq_along(panels)) {
    p <- panels[[i]]
    if (!p$matrix %in% COMPARISON_MATRICES) next
    r_b <- hic_resolution(p$vs, locus_res)
    if (r_b != panel_res[i]) {
      stop("hic.panels[[", p$index, "]] (", p$name, ") at locus ", nomen, ": '", p$hic,
           "' is at ", panel_res[i], " bp but '", p$vs, "' is at ", r_b, " bp. A '",
           p$matrix, "' panel merges the two matrices bin by bin, so both files must ",
           "use the same resolution.", call. = FALSE)
    }
  }

  # Pad the requested region, then snap the window to the COARSEST resolution
  # in use, so that every panel's grid -- whatever its own bin size -- starts
  # and ends on a shared boundary.
  region_size <- reg$end - reg$start
  grid_res <- max(panel_res)
  win_start <- round((reg$start - region_size * locus$padding_fraction) / grid_res) * grid_res
  win_end   <- round((reg$end   + region_size * locus$padding_fraction) / grid_res) * grid_res

  # One crop, shared by every panel and by the track panel below.
  lims <- crop_limits(win_start, win_end, crop_x, crop_y)

  # ---- heatmap panels -------------------------------------------------------
  # Values first, plots second: a `normalize` group needs every member's values
  # before any of their colour scales can be built.
  maps <- lapply(seq_along(panels), function(i) {
    panel_locus_map(panels[[i]], chr_id, win_start, win_end, panel_res[i])
  })
  # Limits are computed from the upper triangle, which is what actually gets drawn.
  own_vals <- lapply(maps, function(m) m[y >= x, value])

  scale_vals <- own_vals
  groups <- unique(na.omit(vapply(panels, function(p) p$normalize, character(1))))
  for (grp in groups) {
    idx <- which(vapply(panels, function(p) identical(p$normalize, grp), logical(1)))
    pooled <- unlist(own_vals[idx], use.names = FALSE)
    for (i in idx) scale_vals[[i]] <- pooled
    message("  normalize group '", grp, "': ",
            paste(vapply(panels[idx], function(p) p$name, character(1)), collapse = ", "))
  }

  plots <- lapply(seq_along(panels), function(i) {
    matriz_horizontal(panels[[i]], maps[[i]], chr_label, panel_res[i], lims, scale_vals[[i]])
  })

  # Only the first panel keeps the x axis and the title; the rest are bare so
  # they stack seamlessly.
  plot_ls <- lapply(seq_along(plots), function(i) {
    if (i == 1) {
      plots[[i]] + theme(
        axis.text.x.top  = element_text(size = 6, margin = margin(b = -1)),
        axis.ticks.x.top = element_line(color = "black"),
        plot.margin = margin(0, 0, 0, 20),
        plot.title  = element_text(size = 6, margin = margin(b = -1)))
    } else {
      plots[[i]] + theme(
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), plot.title = element_blank(),
        plot.margin = margin(0, 0, 0, 0))
    }
  })

  # Every track aligns to the same visible interval as the heatmaps.
  visible_region <- lims$visible
  region_gr <- GRanges(seqnames = chr_label,
                       ranges = IRanges(visible_region$start, visible_region$end))
  window <- visible_region$end - visible_region$start

  # ---- marks ----------------------------------------------------------------
  mk_tracks <- tracks_of("mark")
  bw_signal_dt <- bin_marks(region_gr)
  scale_labels <- data.table()
  if (nrow(bw_signal_dt)) {
    bw_signal_dt[, log_signal := log10(1 + signal)]
    bw_signal_dt[, display_max := max(log_signal, na.rm = TRUE), by = mark]
    # Marks in the same group get one common maximum, so their heights are
    # directly comparable to each other.
    for (group in shared_groups) {
      grp <- unlist(group)
      bw_signal_dt[mark %in% grp, display_max := max(log_signal, na.rm = TRUE)]
    }
    # Attach each mark's computed band and drawing options.
    mk_dt <- rbindlist(lapply(mk_tracks, function(t) data.table(
      mark = t$name, ymin = t$ymin, band = t$height, signal_scale = t$signal_scale,
      crop = t$crop, color = t$color, alpha = t$alpha,
      line_color = t$line_color, line_width = t$line_width, show_scale = t$show_scale)))
    bw_signal_dt <- mk_dt[bw_signal_dt, on = "mark"]

    # signal_scale: auto -> exactly fill the band at this locus
    bw_signal_dt[is.na(signal_scale),
                 signal_scale := band / pmax(display_max, .Machine$double.eps)]
    bw_signal_dt[, ytop := ymin + (log_signal * signal_scale)]
    bw_signal_dt[crop == TRUE, ytop := pmin(ytop, ymin + band)]

    scale_labels <- bw_signal_dt[show_scale == TRUE,
                                 .(max = round(unique(display_max), 2),
                                   ymin = unique(ymin), band = unique(band)), by = mark]
  }

  # ---- structures -----------------------------------------------------------
  st_parts <- lapply(tracks_of("structure"), build_structure_track,
                     chr_label = chr_label,
                     win_start = visible_region$start, win_end = visible_region$end)
  names(st_parts) <- vapply(tracks_of("structure"), function(t) t$name, character(1))

  # ---- annotations ----------------------------------------------------------
  an_tracks <- tracks_of("annotation")
  an_parts <- lapply(an_tracks, build_annotation_track,
                     chr_id = chr_id, region_gr = region_gr, window = window)
  gather <- function(key) {
    dts <- Filter(function(d) !is.null(d) && nrow(d), lapply(an_parts, `[[`, key))
    if (!length(dts)) data.table() else rbindlist(dts, fill = TRUE)
  }
  # Rows of one geometry belonging to one track. Guards the empty case, where
  # the data.table has no columns at all and `track` would not resolve.
  pick <- function(d, nm) if (!nrow(d) || !"track" %in% names(d)) data.table() else d[track == nm]
  an_boxes <- gather("boxes");    an_introns <- gather("introns")
  an_exons <- gather("exons");    an_cds     <- gather("cds")
  an_chev  <- gather("chevrons"); an_tri     <- gather("triangles")
  an_labs  <- gather("labels")

  # =========================================================================
  # TRACK PANEL: layers are added back to front
  # =========================================================================
  tracks_plot <- ggplot()

  # 1) structures, behind everything else
  for (s in tracks_of("structure")) {
    parts <- st_parts[[s$name]]
    # `alpha` and `line` are one value per FEATURE, computed by the builder, and
    # reach the geom untouched through the identity scales added at the end.
    if (nrow(parts$domains)) {
      tracks_plot <- tracks_plot + geom_rect(
        data = parts$domains,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, group = id,
            alpha = alpha, color = line),
        fill = if (identical(s$fill, "none")) NA else s$fill,
        linewidth = s$line_width, linetype = s$line_type, inherit.aes = FALSE)
    }
    if (nrow(parts$arcs)) {
      # filled area and outline are independent, so a track can be either
      if (!identical(s$fill, "none")) {
        tracks_plot <- tracks_plot + geom_polygon(
          data = parts$arcs, aes(x = x, y = y, group = id, alpha = alpha),
          fill = s$fill, inherit.aes = FALSE)
      }
      if (!identical(s$line_color, "none")) {
        tracks_plot <- tracks_plot + geom_path(
          data = parts$arcs, aes(x = x, y = y, group = id, color = line),
          linewidth = s$line_width,
          linetype = s$line_type, inherit.aes = FALSE)
      }
    }
  }

  # 2) bigwig ribbons, one layer per mark so each keeps its own alpha and outline
  if (nrow(bw_signal_dt)) {
    for (t in mk_tracks) {
      d <- bw_signal_dt[mark == t$name]
      if (!nrow(d)) next
      tracks_plot <- tracks_plot + geom_ribbon(
        data = d, aes(x = start, ymin = ymin, ymax = ytop, fill = color, group = mark),
        alpha = t$alpha, color = t$line_color, linewidth = t$line_width)
    }
    # "max" and "0" markers for each mark's own scale
    if (nrow(scale_labels)) {
      tracks_plot <- tracks_plot +
        geom_text(data = scale_labels,
                  aes(x = visible_region$start + window * 0.005,
                      y = ymin + band * 0.88, label = max),
                  inherit.aes = FALSE, size = 2, hjust = 0) +
        geom_text(data = scale_labels,
                  aes(x = visible_region$start + window * 0.005,
                      y = ymin + band * 0.12, label = 0),
                  inherit.aes = FALSE, size = 2, hjust = 0)
    }
  }

  # 3) annotations: intron lines first, then exons, then CDS on top of exons.
  #    Drawn per track so each keeps its own outline color and width.
  for (a in an_tracks) {
    lc <- a$line_color; lw <- a$line_width
    d <- pick(an_introns, a$name)
    if (nrow(d)) tracks_plot <- tracks_plot + geom_segment(
      data = d, aes(x = x, xend = xend, y = y, yend = yend),
      color = lc, linewidth = lw, inherit.aes = FALSE)
    for (part in list(an_exons, an_cds, an_boxes)) {
      d <- pick(part, a$name)
      if (nrow(d)) tracks_plot <- tracks_plot + geom_rect(
        data = d, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = color),
        color = lc, linewidth = lw, inherit.aes = FALSE)
    }
    d <- pick(an_tri, a$name)
    if (nrow(d)) tracks_plot <- tracks_plot + geom_polygon(
      data = d, aes(x = x, y = y, group = id, fill = color),
      color = lc, linewidth = lw, inherit.aes = FALSE)
    d <- pick(an_chev, a$name)
    if (nrow(d)) tracks_plot <- tracks_plot + geom_segment(
      data = d, aes(x = x, xend = xend, y = y, yend = yend),
      color = lc, linewidth = lw, inherit.aes = FALSE)
  }

  # 4) annotation labels, on top of everything
  for (i in seq_along(an_tracks)) {
    a <- an_tracks[[i]]
    d <- pick(an_labs, a$name)
    if (!nrow(d)) next
    lab <- a$labels
    if (isTRUE(lab$repel %||% TRUE)) {
      # The box repulsion may move a label inside, computed once by the builder
      # so the two places cannot drift apart. Pinning it to the band centre
      # would drag every label back onto the features and undo the reserve.
      zone       <- an_parts[[i]]$label_zone
      lab_centre <- zone$centre
      slack      <- zone$slack
      tracks_plot <- tracks_plot + geom_text_repel(
        data = d, aes(x = x, y = y, label = label),
        # nudge_y moves the label while leaving the anchor on the feature, so
        # the leader line keeps pointing at the gene it names.
        nudge_y = d$nudge,
        size = as.numeric(lab$size %||% 2),
        color = as.character(lab$color %||% "black"),
        fontface = as.character(lab$fontface %||% "bold"),
        force = as.numeric(lab$force %||% 15), max.overlaps = Inf,
        # direction: "x" keeps every label at its height and spreads them
        # sideways -- usually what you want above a gene track; "y" stacks
        # them vertically; "both" lets repulsion use either axis.
        direction = as.character(lab$direction %||% "both"),
        min.segment.length = as.numeric(lab$min_segment_length %||% 0),
        segment.size = as.numeric(lab$segment_size %||% 0.2),
        segment.color = as.character(lab$segment_color %||% lab$color %||% "black"),
        box.padding = as.numeric(lab$box_padding %||% 0.2), point.padding = 0.1,
        ylim = c(lab_centre - slack, lab_centre + slack), inherit.aes = FALSE)
    } else {
      # No repulsion: draw the label at the anchor plus its nudge, no leader line.
      tracks_plot <- tracks_plot + geom_text(
        data = d, aes(x = x, y = y + nudge, label = label),
        size = as.numeric(lab$size %||% 2),
        color = as.character(lab$color %||% "black"),
        fontface = as.character(lab$fontface %||% "bold"), inherit.aes = FALSE)
    }
  }

  # 5) scales and theme
  #
  # The visible range is set with coord_cartesian, NOT with scale limits.
  # Scale limits DROP any row that falls outside them, which silently deletes
  # whole features that merely straddle the edge -- a gene wider than the view
  # would lose its intron line entirely while its exons still drew. A coord
  # zoom clips the same geometry instead of discarding it.
  # The identity scales let a layer carry its own colour and alpha per row --
  # that is how one structure track shades every feature by its nesting depth
  # without splitting into one layer per level. Nothing here draws a guide.
  tracks_plot <- tracks_plot +
    scale_fill_identity() +
    scale_color_identity() +
    scale_alpha_identity() +
    scale_x_continuous(expand = c(0, 0),
                       labels = function(x) round(x / (x_axis_cfg$divide %||% 1e6), 2),
                       position = "top") +
    scale_y_continuous(expand = c(0, 0),
                       breaks = track_axis_dt$y, labels = track_axis_dt$label) +
    coord_cartesian(xlim = c(visible_region$start, visible_region$end),
                    ylim = tracks_ylim, expand = FALSE) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_blank(), axis.ticks.y = element_blank(),
      axis.text.y = element_text(size = tracks_label_size, angle = 0, hjust = 1, face = "bold"),
      axis.title.x = element_blank(), axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(0, 0, 0, 1),
      panel.border = element_rect(fill = NA, color = "black"))

  # ---- combine and write ----------------------------------------------------
  final_plot <- wrap_plots(c(plot_ls, list(tracks_plot)), ncol = 1) +
    plot_layout(heights = panel_proportions(panels))

  dir.create(locus$dir, recursive = TRUE, showWarnings = FALSE)
  path.plot <- file.path(locus$dir, paste0(nomen, ".", chr_id, "_", reg$start, "_", reg$end,
                                           ".", grid_res, "._bins.hicplot.", locus$formato))
  if (locus$formato == "png") {
    png(path.plot, width = locus$width, height = locus$height, units = "cm", res = locus$dpi)
  } else {
    pdf(path.plot, width = locus$width / 2.54, height = locus$height / 2.54)
  }
  on.exit(dev.off(), add = TRUE)
  print(final_plot)
  message("  wrote ", path.plot)
  invisible(path.plot)
}

# =============================================================================
# MAIN
# =============================================================================
for (locus in loci) {
  render_locus(locus)
}

message("done in ", round(difftime(Sys.time(), start.time, units = "secs"), 1), " s")
