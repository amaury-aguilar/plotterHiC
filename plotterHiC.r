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

# =============================================================================
# ORDER OVERRIDE
# =============================================================================
# Declaration order IS the layout: the first track written is the top band, the
# first panel the top heatmap, and nothing in the config carries a coordinate.
# That stays true -- `plot.order` only exists so that CHANGING the order does
# not mean moving the blocks. A track entry is thirty lines of settings and
# comments; dragging it three places up to put it under the heatmaps is exactly
# the kind of edit that loses a brace.
#
# So `plot.order` names entries, in the order they should be drawn:
#
#   plot: { order: { tracks: [CTCF, ATAC_WT, TADs], panels: [subtraction] } }
#
# The list may be partial: what it names comes first, in that order, and
# everything else follows in declaration order. Bringing one entry to the top
# is therefore a one-line edit, and the entries not mentioned keep the relative
# order the file already gives them.
#
# A name that belongs to a DISABLED entry is skipped rather than refused --
# turning a track off must never force a second edit somewhere else -- but a
# name that belongs to nothing at all is a typo, and stops the run.
apply_order <- function(items, handles, want, all_handles, what) {
  want <- as.character(unlist(want %||% list()))
  if (!length(want)) return(items)

  dup <- unique(want[duplicated(want)])
  if (length(dup)) {
    stop(what, ": '", paste(dup, collapse = ", "), "' listed more than once",
         call. = FALSE)
  }
  bad <- setdiff(want, all_handles)
  if (length(bad)) {
    stop(what, ": '", paste(bad, collapse = ", "), "' is not defined. Use one of: ",
         paste(all_handles, collapse = ", "), call. = FALSE)
  }
  amb <- intersect(want, handles[duplicated(handles)])
  if (length(amb)) {
    stop(what, ": '", paste(amb, collapse = ", "), "' matches more than one entry;",
         " give them distinct names before ordering by name", call. = FALSE)
  }

  idx <- match(want, handles)          # NA = declared but disabled: nothing to place
  idx <- idx[!is.na(idx)]
  items[c(idx, setdiff(seq_along(items), idx))]
}

plot_order <- cfg$plot$order %||% list()
if (length(plot_order)) {
  bad <- setdiff(names(plot_order), c("tracks", "panels"))
  if (length(bad)) {
    stop("plot.order: '", paste(bad, collapse = ", "), "' is not valid. Use one of: ",
         "tracks, panels", call. = FALSE)
  }
}

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
# What a panel answers to in `plot.order.panels`: its own `name`, or -- for the
# panel left unnamed because its title is drawn by its legend -- its matrix
# type, so `panels: [subtraction]` reaches it without inventing a name for it.
panel_handles <- vapply(seq_along(panels_cfg), function(i) {
  ef <- merge_settings(panels_cfg[[i]], panel_defaults$flat, panel_presets,
                       paste0("hic.panels[[", i, "]]"))
  nm <- as.character(ef$name %||% "")
  if (nzchar(nm)) nm else as.character(ef$matrix %||% "log2")
}, character(1))

panels_cfg <- panels_cfg[panel_on]
if (!length(panels_cfg)) {
  stop("every entry of hic.panels has 'enabled: false' in ", cfg_path,
       ": at least one panel must be on", call. = FALSE)
}
panels_cfg <- apply_order(panels_cfg, panel_handles[panel_on], plot_order$panels,
                          panel_handles, "plot.order.panels")

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
track_names_all <- names(track_list)     # for plot.order, disabled ones included
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

    # n_bins and bin_size answer the same question two ways (see mark_n_bins),
    # so one track cannot use both -- but the test is on what THIS block
    # writes, not on the flattened settings: a `tracks.defaults` that sets
    # n_bins for every mark must not stop one track from asking for a bin_size.
    # Whichever the entry itself names wins, and cancels the other's default.
    own_nb <- has_key(raw, "n_bins")
    own_bs <- has_key(raw, "bin_size")
    if (own_nb && own_bs) {
      stop("tracks.list.", tn, ": set either n_bins (a fixed count of bins) or ",
           "bin_size (a fixed bin width in bp), not both", call. = FALSE)
    }
    c(base, list(
      color      = as.character(opt(e, d, "color", "grey50")),
      alpha      = as.numeric(opt(e, d, "alpha", 0.9)),
      line_color = as.character(opt(e, d, "line_color", "black")),
      line_width = as.numeric(opt(e, d, "line_width", 0.2)),
      signal_scale = ss,
      # HOW FINELY THE BIGWIG IS BINNED -- three ways to say it, and at most
      # one of them per track. See mark_n_bins().
      #   n_bins: 400     a fixed COUNT across the window, whatever its width
      #   bin_size: 200   a fixed WIDTH in bp, however many that takes
      #   n_bins: auto    one bin per output pixel (the default)
      n_bins     = { nb <- if (own_bs) "auto" else opt(e, d, "n_bins", "auto")
                     if (is.character(nb) && tolower(nb) == "auto") NA_real_
                     else as.numeric(nb) },
      bin_size   = { bs <- if (own_nb) NA_real_ else num_or_na(opt(e, d, "bin_size", NA))
                     if (!is.na(bs) && bs <= 0) {
                       stop("tracks.list.", tn, ": bin_size must be a positive ",
                            "number of bp, got ", bs, call. = FALSE)
                     }
                     bs },
      # how the bases inside one bin are reduced to the one value drawn --
      # a genome browser's "windowing function". See bin_marks().
      windowing  = { wf <- tolower(as.character(opt(e, d, "windowing", "mean")))
                     if (!wf %in% c("mean", "max")) {
                       stop("tracks.list.", tn, ": windowing must be 'mean' or 'max', got '",
                            wf, "'", call. = FALSE)
                     }
                     wf },
      crop       = as.logical(opt(e, d, "crop", TRUE)),
      show_scale = as.logical(opt(e, d, "show_scale", TRUE))))

  } else if (ty == "structure") {
    # height_span: the anchor span (bp) whose arc exactly fills the band.
    # "auto" (the default) uses the widest structure in view instead.
    hs <- opt(e, d, "height_span", "auto")
    hs <- if (is.character(hs) && tolower(hs) == "auto") 0 else as.numeric(hs)
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
      # HOW TALL AN ARC IS (style: arc). An arc is drawn as a half ELLIPSE, so
      # its width and its height are set independently: the feet always sit on
      # the anchors, and the apex is whatever this mapping asks for.
      #
      #   span  height is proportional to the loop's length. Faithful, and
      #         unusable as soon as one locus mixes 20 kb loops with 5 Mb ones:
      #         the long one takes the band and the local ones are a hairline.
      #   log   height is proportional to log(length), so a 200x range of loop
      #         lengths still spreads over the band. This is what keeps the
      #         local loops readable next to a long-range one.
      #   rank  height by ORDER of length, evenly spread over the band: every
      #         arc is distinguishable whatever the lengths are, at the cost of
      #         the height no longer meaning a length. Pair it with height_min.
      height_by = as.character(opt(e, d, "height_by", "span")),
      # The span drawn one full band tall -- in bp (height_span) or as a
      # fraction of the visible window (height_span_frac), the larger winning.
      # auto/0 = the widest structure in view, so nothing exceeds the band.
      height_span      = hs,
      height_span_frac = as.numeric(opt(e, d, "height_span_frac", 0)),
      # log only: the span at the BOTTOM of the scale. auto/0 = the shortest in
      # view. Raise it to stop one unusually short loop from setting the floor.
      height_span_min  = as.numeric(opt(e, d, "height_span_min", 0)),
      # The mapping is rescaled into [height_min, height_max], both fractions of
      # the band: the floor keeps the shortest arc off the baseline, and the
      # ceiling keeps the tallest one from touching the divider drawn at the top
      # of the track, where it reads as part of the frame rather than as a loop.
      height_min = as.numeric(opt(e, d, "height_min", 0)),
      height_max = as.numeric(opt(e, d, "height_max", 1)),
      # Safety net for `height_by: span` with an explicit height_span, the one
      # combination that can ask for an arc taller than its band: TRUE cuts it
      # flat at the top, FALSE lets it draw over the tracks above. It does
      # nothing under log/rank, where the mapping cannot exceed the band --
      # and cutting is ugly by nature: two cut arcs share one ceiling and read
      # as a single horizontal line notched wherever they meet.
      crop = as.logical(opt(e, d, "crop", TRUE)),
      # Keep only structures whose anchor span is in [min_span, max_span], in bp
      # or as a fraction of the visible window; the larger of each pair wins,
      # 0 = off. The _frac form is what travels between loci of different sizes:
      # a loop many times longer than the window has its apex flat across the
      # whole panel and no anchor to show for it, so it reads as a horizontal
      # rule over the track however its height was computed. That is a job for
      # a filter, not for the height mapping.
      min_span  = as.numeric(opt(e, d, "min_span", 0)),
      max_span  = as.numeric(opt(e, d, "max_span", 0)),
      min_span_frac = as.numeric(opt(e, d, "min_span_frac", 0)),
      max_span_frac = as.numeric(opt(e, d, "max_span_frac", 0)),
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
      # NON-STANDARD BED: which columns to read, and what they mean. Naming
      # them bypasses rtracklayer, which reads a BED by position and meaning and
      # chokes on a file that keeps its own columns past the sixth.
      #   select:  [1, 2, 3, 4, 5, 6]
      #   columns: [chr, start, end, name, score, strand]
      select  = opt(e, d, "select"),
      columns = opt(e, d, "columns"),
      # NA = take the color from the BED's itemRgb column
      color      = as.character(opt(e, d, "color", NA_character_)),
      # ...or drop both and colour each feature by a NUMBER instead:
      #   none   flat `color`, else the BED's itemRgb, else black
      #   score  the BED's score column (5) through the diverging ramp in
      #          `color_scale`. A DEG bed stores a signed log2FC there, so the
      #          track then shows direction and magnitude rather than class.
      # `score` overrides `color` and itemRgb -- a flat colour and a gradient
      # cannot both be honoured, and the one the config asked for wins.
      color_by    = { cb <- tolower(as.character(opt(e, d, "color_by", "none")))
                      if (!cb %in% c("none", "score")) {
                        stop("tracks.list.", tn, ": 'color_by: ", cb, "' is not",
                             " valid. Use one of: none, score", call. = FALSE)
                      }
                      cb },
      color_scale = opt_block(e, d, "color_scale"),
      line_color = as.character(opt(e, d, "line_color", "black")),
      line_width = as.numeric(opt(e, d, "line_width", 0.2)),
      utr_thinner = as.logical(opt(e, d, "utr_thinner", TRUE)),
      utr_frac    = as.numeric(opt(e, d, "utr_frac", 0.45)),
      # grow features narrower than this fraction of the window so single-bp
      # exons stay visible; position and overall span are preserved
      min_exon_frac = as.numeric(opt(e, d, "min_exon_frac", 0)),
      # FLOOR ON THE WHOLE FEATURE: a gene narrower than this is drawn zoomed
      # about its own midpoint until it reaches the floor, exons and all. In bp
      # (min_feature_width) or as a fraction of the visible window
      # (min_feature_frac); the larger of the two wins. 0 = off.
      # Unlike min_exon_frac, this DOES stretch the drawn span of the gene: a
      # 2 kb gene in a 1.5 Mb view is a hairline otherwise.
      min_feature_frac  = as.numeric(opt(e, d, "min_feature_frac", 0)),
      min_feature_width = as.numeric(opt(e, d, "min_feature_width", 0)),
      # chevrons = repeated ">" marks along the body; triangle = one arrowhead
      direction  = as.character(opt(e, d, "direction", "triangle")),
      chevron_n  = as.numeric(opt(e, d, "chevron_n", 6)),
      # Length of the arrowhead, as a fraction of the VISIBLE WINDOW -- so it is
      # the same size whatever the gene's length, and a short gene still gets a
      # full-sized head. Raise it to make small features carry the track.
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
# The one place the track order is decided, so the bands below and everything
# that reads them follow from it.
tracks <- apply_order(tracks, names(tracks), plot_order$tracks,
                      track_names_all, "plot.order.tracks")

# --- compute each band from the resolved order -------------------------------
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
    # Both sides are whole-chromosome tables, and the join below keeps only what
    # falls on the grid -- so trim to the window BEFORE the outer merge rather
    # than merging millions of bin pairs and throwing nearly all of them away.
    # The cached tables themselves are left untouched, other loci still need them.
    in_win <- function(d) d[x >= win_start & x <= win_end & y >= win_start & y <= win_end]
    m <- merge(in_win(a$dt), in_win(b$dt), by = c("x", "y"), all = TRUE,
               suffixes = c("", "_vs"))
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

  # ...but only INSIDE the chromosome. A grid anchor past either end is not an
  # empty bin, it is not a bin at all: `padding_fraction` adds context that does
  # not exist at a telomere, and a whole-chromosome locus is usually written
  # past the end on purpose. Filling those with the neutral value paints them
  # the colour of "average" -- which on a log2 ramp is a solid orange field
  # wrapped around the data triangle, and reads as data. Left NA, the panel's
  # own `na_color` decides what the outside of the chromosome looks like.
  size <- chr.list[chr == chr_id, size]
  if (length(size)) m[x < 0 | y < 0 | x >= size | y >= size, value := NA_real_]
  m[, .(x, y, value)]
}

# =============================================================================
# HEATMAP RENDERING
# =============================================================================

# The four corners of EVERY matrix bin, rotated -45 degrees so the diagonal of
# the contact matrix becomes the horizontal genomic axis.
#
# Vectorised over the whole table on purpose. Called one bin at a time -- which
# is what `by = .(x, y)` did -- this costs about 120 us per bin in data.table
# group overhead alone, so a 5 Mb window at 10 kb (350k bins in the upper
# triangle) spent 43 s per panel building corners for a map that renders in 1 s.
to_polygon <- function(x, y, counts, resol) {
  half <- resol / 2
  n    <- length(x)
  # Corner order within each bin is 1..4 around the square, as before, so
  # geom_polygon still closes it: (-,-) (+,-) (+,+) (-,+).
  cx <- rep(x, each = 4L) + c(-half,  half, half, -half)
  cy <- rep(y, each = 4L) + c(-half, -half, half,  half)
  data.table(u = (cx + cy) / sqrt(2), v = (cy - cx) / sqrt(2),
             counts = rep(counts, each = 4L),
             id     = rep(seq_len(n), each = 4L),
             corner = rep_len(1:4, 4L * n))
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
  bins <- plot_df[y >= x]

  # ...and of that triangle, only what the crop keeps. A bin whose rotated
  # square lies entirely outside the axis limits is dropped by ggplot anyway,
  # at the scale, after its corners have been built and carried through the
  # whole pipeline -- and with the usual crop that is half the triangle. The
  # test is on the bin's own extent, so a square straddling a limit is kept and
  # the picture is unchanged. Colour limits are NOT affected: `scale_vals` is
  # computed by the caller from the full triangle.
  s2 <- sqrt(2)
  bins <- bins[(x + y + resol) / s2 >= lims$u_min &
               (x + y - resol) / s2 <= lims$u_max &
               (y - x - resol) / s2 <= lims$v_max]

  polygonos <- to_polygon(bins$x, bins$y, bins$counts, resol)

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

# How many bins a mark gets across the visible region -- the single most
# important number for how a sharp mark reads, and the one no scaling knob
# downstream can compensate for.
#
# A CTCF site is a few hundred bp. Averaged into a 12.5 kb bin it is diluted
# ~40x into the local background, and on chr15:55-60 Mb that takes the
# summit-to-background ratio of the CN sample from 6.6 (at 1 kb bins) down to
# 1.7. At 1.7 nothing downstream can help: `signal_scale: auto` divides the
# band by the maximum, so a background sitting at 1/1.7 of that maximum fills
# 59% of the band and the track reads as a noisy slab. The dynamic range is
# gone here, in the binning, before signal_scale ever sees it.
#
# THE THREE MODES. A track picks one:
#
#   n_bins: 400     A fixed COUNT. The window is cut into this many bins
#                   whatever its width, so the bin is 12.5 kb at 5 Mb and
#                   250 kb at 100 Mb -- the resolution changes under you as
#                   soon as you change the locus, which is what makes one
#                   signal_scale that worked at one zoom wrong at the next.
#
#   bin_size: 200   A fixed WIDTH in bp, and as many bins as that takes. The
#                   resolution is now a property of the track rather than of
#                   the window: 200 bp resolves single CTCF sites at any zoom,
#                   and the cost is the bin COUNT, which grows with the window
#                   (a whole chromosome at 200 bp is half a million bins).
#
#   n_bins: auto    One bin per pixel of the finished figure (the default).
#                   This is what a genome browser does -- IGV re-bins to the
#                   window it is drawing, which is the reason a track looks
#                   sharper there than under a fixed n_bins -- and it is the
#                   finest binning that can actually be seen: past one bin per
#                   pixel the extra bins land on top of each other.
#
# ...which is also the ceiling, and `bin_size` is the one mode that can ask to
# go past it: 200 bp is 25k bins over 5 Mb and 390k over a chromosome. Those
# extra bins are not just invisible, they are ruinous -- the whole cost lands
# in the PNG rasterizer drawing a ribbon outline of three quarters of a million
# vertices, and the same figure that takes 14 s at one bin per pixel had not
# finished after 13 minutes.
#
# Capping is exact rather than an approximation, because both windowing
# functions survive the collapse when the bins are equal width: the mean of
# equal-width means is the overall mean, and the max of maxes is the overall
# max. So a capped `bin_size` draws precisely what `n_bins: <cap>` would, and
# the only thing lost is resolution nothing could have displayed.
mark_n_bins <- function(t, locus, region_gr) {
  px <- if (identical(locus$formato, "png")) locus$width / 2.54 * locus$dpi else 2000
  px <- max(200, min(round(px), 20000))

  want <- if (!is.na(t$bin_size)) max(1, round(width(region_gr) / t$bin_size))
          else if (!is.na(t$n_bins)) max(1, round(t$n_bins))
          else return(px)
  if (want <= px) return(want)

  # Say so rather than quietly redrawing the track at a resolution the config
  # did not ask for. What is worth flagging is the bin WIDTH that comes back,
  # because signal_scale is tuned against a bin width: a number that filled the
  # band at 200 bp is wrong at 42 kb.
  #
  # Deliberately no advice to reach for `windowing: max` here. It rescues a bin
  # a few times wider than the peaks, but this cap lands on bins hundreds of
  # times wider, where essentially every bin contains a summit -- so max is
  # high everywhere, and the track saturates into a slab instead of the flat
  # one `mean` gives. Neither is a picture of peaks; the window is simply too
  # wide to have any.
  message("  ", t$name, ": ", if (!is.na(t$bin_size)) paste0("bin_size ", t$bin_size, " bp")
          else paste0("n_bins ", want),
          " over ", round(width(region_gr) / 1e6, 1), " Mb is ", want,
          " bins, more than this figure has pixels -- drawing ", px,
          " of ", round(width(region_gr) / px), " bp. Re-tune signal_scale for ",
          "this width, or set it to 'auto'.")
  px
}

# One value per bin, one row per bin.
#
# `windowing` decides how the bases inside a bin become that one value, and for
# a sharp mark it matters more than any scaling knob downstream.
#
#   mean  per-base average. Weighted by the bases each interval contributes,
#         which is what "mean" has to mean for a bigwig: these files are
#         run-length encoded, so one 400 bp interval over a peak and one 20 kb
#         interval of background are a single row each. An unweighted mean over
#         overlapping intervals -- what this used to do -- counts those two
#         equally and reports a value the region never has.
#
#   max   the tallest value in the bin, i.e. the summit survives whatever the
#         bin size. This is the one to reach for on a low signal-to-noise
#         sample: a CTCF site is a few hundred bp, so at a 2 kb bin a mean
#         spreads a summit over 5-10x its own width and buries it in the local
#         background. On chr15:55-60 Mb that is the whole difference for the CN
#         ChIP -- peak-to-background 3.8 under `mean` against 14.4 under `max`,
#         which is a background slab filling 30% of the band against 14%.
#         `mean` stays the default: it is the honest summary of a broad mark,
#         and `max` on a noisy sample also lifts the noise spikes.
bin_marks <- function(region_gr, locus) {
  mk_tracks <- tracks_of("mark")
  if (!length(mk_tracks)) return(data.table())
  rbindlist(lapply(mk_tracks, function(t) {
    bins <- tile(region_gr, n = mark_n_bins(t, locus, region_gr))[[1]]
    bw   <- import(t$file, which = region_gr)
    sig  <- numeric(length(bins))
    if (length(bw)) {
      ov <- findOverlaps(bins, bw)
      if (length(ov)) {
        qh <- queryHits(ov); sh <- subjectHits(ov)
        if (t$windowing == "max") {
          agg <- tapply(bw$score[sh], qh, max)
          sig[as.integer(names(agg))] <- agg
        } else {
          # bases of each interval that actually fall inside the bin
          wt  <- as.numeric(width(pintersect(bins[qh], bw[sh])))
          agg <- tapply(bw$score[sh] * wt, qh, sum)
          sig[as.integer(names(agg))] <- agg
          # ...divided by the WHOLE bin, not just its covered part: a gap in the
          # bigwig is zero coverage, not missing data, and must pull the mean down.
          sig <- sig / as.numeric(width(bins))
        }
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

  win <- win_end - win_start
  lo_span <- max(s$min_span, s$min_span_frac * win)
  hi_span <- max(s$max_span, s$max_span_frac * win)
  if (lo_span > 0) sc <- sc[(end - start) >= lo_span]
  if (hi_span > 0) sc <- sc[(end - start) <= hi_span]
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
  # band's baseline.
  #
  # Deepest last, so when the arcs are filled the small nested ones are not
  # buried under the wide one that contains them.
  ord <- order(sc$depth)
  sc <- sc[ord]; alp <- alp[ord]
  if (length(line) > 1) line <- line[ord]
  res <- max(3, s$arc_resolution)

  # HOW TALL EACH ARC IS DRAWN. A semicircle ties the apex to the span, which is
  # the whole problem: put a 5 Mb loop and a 20 kb one in the same window and
  # the short one is 0.4% of the band. So the height is decided HERE, from the
  # span, and the shape below is stretched to match it -- the arc is a half
  # ellipse, its feet on the anchors and its apex wherever this says.
  spans <- sc$end - sc$start
  # The span that reaches the top of the band. auto = the widest in view, so
  # nothing can exceed the band whatever the locus holds.
  ref <- max(s$height_span, s$height_span_frac * (win_end - win_start))
  if (ref <= 0) ref <- max(spans)
  frac <- switch(s$height_by,
    # Proportional. Honest, and only readable while the loops of one locus are
    # within a few fold of each other.
    span = spans / ref,
    # Proportional to the LOG of the span. A locus holding 20 kb, 300 kb and
    # 5 Mb loops spreads them over the band instead of pinning the first two to
    # the baseline; the ordering is still by length. The bottom of the scale is
    # the shortest arc in view unless height_span_min names it.
    log = {
      lo <- if (s$height_span_min > 0) s$height_span_min else min(spans)
      if (ref <= lo) rep(1, length(spans))
      else pmin(1, pmax(0, log(pmax(spans, 1) / lo) / log(ref / lo)))
    },
    # Order only, evenly spread. Immune to the span distribution: n distinct
    # lengths get n evenly spaced heights, so nothing is ever squashed. Ties
    # share a height, and the shortest sits on the floor -- give it a
    # height_min or it lies flat on the baseline.
    rank = {
      r <- frank(spans, ties.method = "dense")
      if (max(r) <= 1) rep(1, length(r)) else (r - 1) / (max(r) - 1)
    },
    stop("tracks.list.", s$name, ": 'height_by: ", s$height_by,
         "' is not valid. Use one of: span, log, rank", call. = FALSE))

  # The mapping is RESCALED into [height_min, height_max] rather than clamped,
  # so raising the floor lifts the smallest arcs clear of the baseline without
  # collapsing the differences between them, and lowering the ceiling gives the
  # longest ones headroom under the track divider.
  hmin <- max(0, min(1, s$height_min))
  hmax <- max(hmin, min(1, s$height_max))
  peak <- (hmin + (hmax - hmin) * frac) * s$height

  out <- lapply(seq_len(nrow(sc)), function(j) {
    x1 <- sc$start[j]; x2 <- sc$end[j]
    if (x2 <= x1) return(NULL)
    xs  <- seq(x1, x2, length.out = res)
    rad <- (x2 - x1) / 2; xc <- (x1 + x2) / 2
    # Half ellipse: same shape as the old semicircle, but scaled vertically to
    # the height decided above instead of to the radius.
    y <- s$ymin + peak[j] * sqrt(pmax(0, 1 - ((xs - xc) / rad)^2))
    if (s$crop && peak[j] > s$height) {
      # Only `height_by: span` can get here. Insert the two x where the ellipse
      # crosses the top of the band, so the cut starts exactly there instead of
      # at the nearest sampled point.
      dx <- rad * sqrt(pmax(0, 1 - (s$height / peak[j])^2))
      xs <- sort(c(xs, xc - dx, xc + dx))
      y  <- pmin(s$ymin + peak[j] * sqrt(pmax(0, 1 - ((xs - xc) / rad)^2)), s$ymax)
    }
    data.table(x = xs, y = y,
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

# Per-feature colour from a signed number -- a log2FC, a delta, a correlation.
# A diverging ramp with `mid` pinned to `midpoint`, so the neutral value keeps
# the neutral colour whatever happens to be in view.
#
# The extent is symmetric about the midpoint by construction (the wider limit
# sets it): a ramp whose neutral colour sits off centre reads as a sequential
# one and quietly shows a mild change as a strong one. Values past the limits
# are squished onto the ends rather than dropped, so `limits` is what decides
# how much of the range the eye actually resolves -- set it wide and everything
# real is pale, set it tight and every gene saturates.
#
# The ends of the ramp are resolved once, here, because the legend has to draw
# the same scale the features were coloured with -- `limits: auto` reads the
# data in view, so the two would drift apart if each worked it out on its own.
score_limits <- function(v, cs, track_name) {
  midpoint <- as.numeric(cs$midpoint %||% 0)
  lim <- cs$limits %||% "auto"
  if (is.character(lim) && tolower(lim[1]) == "auto") {
    # Symmetric around the midpoint, from what is in view. Convenient, and the
    # reason to avoid it in a series: the same log2FC then draws a different
    # colour in every figure.
    r <- if (all(is.na(v))) NA_real_ else max(abs(v - midpoint), na.rm = TRUE)
    if (!is.finite(r) || r == 0) r <- 1
    return(c(midpoint - r, midpoint + r))
  }
  lim <- as.numeric(unlist(lim))
  if (length(lim) != 2 || anyNA(lim)) {
    stop("tracks.list.", track_name, ": 'color_scale.limits' must be two",
         " numbers or 'auto'", call. = FALSE)
  }
  lim
}

score_colors <- function(v, cs, track_name, lim = score_limits(v, cs, track_name)) {
  low    <- as.character(cs$low  %||% "blue3")
  mid    <- as.character(cs$mid  %||% "white")
  high   <- as.character(cs$high %||% "red3")
  na_col <- as.character(cs$na_color %||% "grey70")
  midpoint <- as.numeric(cs$midpoint %||% 0)

  if (all(is.na(v))) {
    warning("tracks.list.", track_name, ": 'color_by: score' but this BED has no",
            " usable score column; every feature takes color_scale.na_color",
            call. = FALSE)
    return(rep(na_col, length(v)))
  }

  t <- scales::rescale_mid(v, to = c(0, 1), from = range(lim), mid = midpoint)
  out <- scales::div_gradient_pal(low, mid, high, space = "Lab")(pmin(1, pmax(0, t)))
  out[is.na(out)] <- na_col
  out
}

# The colour bar for a `color_by: score` track: a strip of the same ramp, placed
# beside the track panel the way the heatmaps place theirs. It is a plot of its
# own rather than a ggplot guide because the track panel spends its fill scale
# on `scale_fill_identity()` -- that is what lets one layer carry a colour per
# feature -- and a plot has exactly one scale per aesthetic.
#
# Returns an inset, or NULL when the track is not on a gradient, drew nothing
# here, or turned the legend off.
build_score_legend <- function(a, lim, ylim) {
  if (!identical(a$color_by, "score") || is.null(lim)) return(NULL)
  cs <- a$color_scale
  lg <- cs$legend %||% list()
  if (!isTRUE(lg$show %||% TRUE)) return(NULL)

  # Strips rather than an axis: the tick labels are drawn as geometry inside the
  # bar's own coordinates, so they take a FIXED share of the inset. Left to a
  # ggplot axis, the text width is whatever the numbers happen to be, and it
  # eats the inset from the right until the ramp itself is a hairline.
  edges  <- seq(lim[1], lim[2], length.out = 129)
  strips <- data.table(ymin = utils::head(edges, -1), ymax = utils::tail(edges, -1))
  strips[, col := score_colors((ymin + ymax) / 2, cs, a$name, lim)]
  # A bar this small is a few strips per pixel, and abutting rects leave
  # antialiased seams that read as banding. Overlapping them by half a strip
  # covers the seam; the coord clips whatever overflows the ends.
  strips[, ymax := ymax + (edges[2] - edges[1]) / 2]

  brks <- if (!is.null(lg$breaks)) as.numeric(unlist(lg$breaks))
          else seq(lim[1], lim[2], length.out = max(2, as.numeric(lg$n_breaks %||% 3)))
  brks <- brks[brks >= min(lim) & brks <= max(lim)]
  digits <- as.numeric(lg$digits %||% 2)
  bw <- 0.42                                  # ramp width, the rest is labels

  bar <- ggplot() +
    geom_rect(data = strips,
              aes(xmin = 0, xmax = bw, ymin = ymin, ymax = ymax, fill = col)) +
    scale_fill_identity()
  if (isTRUE(lg$ticks %||% TRUE)) {
    bar <- bar + geom_segment(data = data.table(y = brks),
                              aes(x = bw, xend = bw * 1.2, y = y, yend = y),
                              linewidth = 0.2, colour = "black")
  }
  # Sizes stay in pt as everywhere else in the config; geom_text wants mm.
  bar <- bar +
    geom_text(data = data.table(y = brks, lab = as.character(signif(brks, digits))),
              aes(x = bw * 1.45, y = y, label = lab), hjust = 0,
              size = as.numeric(lg$text_size %||% 3.5) / .pt) +
    coord_cartesian(xlim = c(0, 1), ylim = range(lim), expand = FALSE, clip = "off") +
    labs(title = as.character(lg$title %||% "score")) +
    theme_void() +
    theme(plot.title = element_text(size = as.numeric(lg$title_size %||% 4),
                                    hjust = 0, face = "bold", margin = margin(b = 1)),
          plot.margin = margin(0, 0, 0, 0))

  # Placement is in fractions of the TRACK PANEL, and x is normally negative:
  # the bar lives in the margin to its left, on the same column as the heatmap
  # legends. `auto` on y centres it on the band of its own track, which is what
  # ties the scale to the features it explains without hardcoding a number that
  # breaks the moment a track is enabled above it.
  pos <- lg$position %||% list(-0.09, "auto")
  x   <- as.numeric(pos[[1]])
  yp  <- pos[[2]]
  ycen <- if (is.character(yp) && tolower(yp) == "auto") {
            (a$ycen - ylim[1]) / diff(ylim)
          } else as.numeric(yp)
  w <- as.numeric(lg$width %||% 0.018)
  h <- if (is.null(lg$height) || identical(tolower(as.character(lg$height)), "auto")) {
         a$height / diff(ylim)          # as tall as the track's own band
       } else as.numeric(lg$height)

  inset_element(bar, left = x - w / 2, bottom = ycen - h / 2,
                right = x + w / 2, top = ycen + h / 2,
                align_to = "panel", clip = FALSE, on_top = TRUE)
}

# The BED columns this reader understands, in their standard order.
BED_FIELDS <- c("chr", "start", "end", "name", "score", "strand", "itemRgb")

# A plain BED reader: takes the columns it is told to take and asks nothing of
# the rest. No blocks, so no exon structure -- `style: model` falls back to
# boxes on its own.
read_bed_plain <- function(file, sel, cols, track_name) {
  dt <- fread(file, sep = "\t", header = FALSE, fill = TRUE, quote = "",
              showProgress = FALSE, colClasses = "character")
  # A '#chrom ...' header, and any track/browser line, are ordinary data rows
  # to fread. Dropping them here is what lets the file keep its header.
  dt <- dt[!grepl("^(#|track|browser)", dt[[1]])]
  if (!nrow(dt)) {
    stop("tracks.list.", track_name, ": no data rows in ", file, call. = FALSE)
  }

  if (is.null(sel) && is.null(cols)) {
    # The standard meaning of the first columns, which is all a BED guarantees.
    sel  <- seq_len(min(ncol(dt), 6L))
    cols <- BED_FIELDS[sel]
  } else if (is.null(sel)) {
    sel <- seq_along(cols)
  } else if (is.null(cols)) {
    cols <- BED_FIELDS[seq_along(sel)]
  }
  if (length(sel) != length(cols)) {
    stop("tracks.list.", track_name, ": `select` and `columns` must name the",
         " same number of fields (", length(sel), " vs ", length(cols), ")",
         call. = FALSE)
  }
  bad <- setdiff(cols, BED_FIELDS)
  if (length(bad)) {
    stop("tracks.list.", track_name, ": 'columns: ", paste(bad, collapse = ", "),
         "' is not valid. Use one of: ", paste(BED_FIELDS, collapse = ", "),
         call. = FALSE)
  }
  if (!all(c("chr", "start", "end") %in% cols)) {
    stop("tracks.list.", track_name, ": `columns` must include chr, start and end",
         call. = FALSE)
  }
  if (max(sel) > ncol(dt)) {
    stop("tracks.list.", track_name, ": `select` asks for column ", max(sel),
         " but ", file, " has ", ncol(dt), call. = FALSE)
  }

  take <- function(nm) if (nm %in% cols) dt[[sel[match(nm, cols)]]] else NULL
  st <- suppressWarnings(as.numeric(take("start")))
  en <- suppressWarnings(as.numeric(take("end")))
  keep <- !is.na(st) & !is.na(en)
  if (!any(keep)) {
    stop("tracks.list.", track_name, ": no row of ", file,
         " has numeric start/end in the selected columns", call. = FALSE)
  }

  strd <- take("strand")
  strd <- if (is.null(strd)) rep("*", sum(keep)) else {
    s <- strd[keep]; s[!s %in% c("+", "-")] <- "*"; s
  }
  # BED is 0-based half-open, the same shift rtracklayer applies, so a feature
  # keeps the coordinates it would have had through the importer.
  gr <- GRanges(seqnames = as.character(take("chr"))[keep],
                ranges  = IRanges(start = st[keep] + 1, end = en[keep]),
                strand  = strd)
  if ("name"  %in% cols) mcols(gr)$name  <- as.character(take("name"))[keep]
  if ("score" %in% cols) mcols(gr)$score <- suppressWarnings(as.numeric(take("score")))[keep]
  if ("itemRgb" %in% cols) mcols(gr)$itemRgb <- as.character(take("itemRgb"))[keep]
  gr
}

# Read one annotation track's file.
#
# rtracklayer's importer is the right tool for a real BED -- it is the only
# thing here that understands BED12 blocks, which `style: model` is built on --
# but it reads columns by POSITION AND MEANING: column 7 must be thickStart, an
# integer. A table that carries its own extra columns (a GWAS bed with the trait
# in column 7) therefore does not fail as a malformed BED; it dies inside scan()
# with "expected 'an integer', got '<some word>'". So the plain reader above
# stands behind it.
#
# Naming `select` (and optionally `columns`) skips the importer altogether,
# which is the honest thing to do for a file that was never a BED past column 6:
#   select:  [1, 2, 3, 4, 5, 6]
#   columns: [chr, start, end, name, score, strand]
load_annotation <- function(a, region_gr) {
  sel  <- if (!is.null(a$select))  as.integer(unlist(a$select))    else NULL
  cols <- if (!is.null(a$columns)) as.character(unlist(a$columns)) else NULL
  if (!is.null(sel) || !is.null(cols)) {
    return(read_bed_plain(a$file, sel, cols, a$name))
  }
  gr <- tryCatch(import(a$file, which = region_gr), error = function(e) {
    warning("tracks.list.", a$name, ": rtracklayer could not read this file as a",
            " BED (", conditionMessage(e), "). Reading the first 6 columns as",
            " chr/start/end/name/score/strand instead -- set `select:` and",
            " `columns:` if that is not what they hold. BED12 blocks are lost on",
            " this path, so `style: model` will draw boxes.", call. = FALSE)
    NULL
  })
  # The importer subsets by region itself; the plain reader reads the file and
  # lets the caller's chromosome and overlap filters do it.
  if (is.null(gr)) read_bed_plain(a$file, NULL, NULL, a$name) else gr
}

# Build every geometry a single annotation track needs.
build_annotation_track <- function(a, chr_id, region_gr, window) {
  empty <- data.table()
  gr <- load_annotation(a, region_gr)
  gr <- gr[seqnames(gr) == paste0("chr", chr_id)]
  gr <- subsetByOverlaps(gr, region_gr)
  if (!length(gr)) {
    warning("tracks.list.", a$name, ": no features in the visible region", call. = FALSE)
    return(list(boxes = empty, introns = empty, exons = empty, cds = empty,
                chevrons = empty, triangles = empty, labels = empty,
                label_zones = data.table(zone = 1L, centre = a$ycen,
                                         slack = a$height / 2)))
  }

  mc <- mcols(gr)
  nms  <- if (!is.null(mc$name))  as.character(mc$name) else rep(NA_character_, length(gr))
  scr  <- if (!is.null(mc$score)) as.numeric(mc$score)  else rep(NA_real_, length(gr))

  # Per-feature colour, in order of precedence: the gradient over the score
  # column if one was asked for, then an explicit track colour, then the BED's
  # own itemRgb, then black. Every glyph the track draws -- exon, CDS, box,
  # arrowhead -- reads this one vector, so the choice is made once here.
  # Kept for the legend, which has to draw the ramp these colours came from.
  score_lims <- if (identical(a$color_by, "score")) {
                  score_limits(scr, a$color_scale, a$name)
                } else NULL
  col <- if (identical(a$color_by, "score")) {
           score_colors(scr, a$color_scale, a$name, score_lims)
         } else if (!is.na(a$color)) {
           rep(a$color, length(gr))
         } else if (!is.null(mc$itemRgb)) {
           vapply(mc$itemRgb, rgb_from_string, character(1))
         } else rep("#000000", length(gr))
  strd <- as.character(strand(gr))

  n_feat <- length(gr)
  tstart <- start(gr); tend <- end(gr)          # TRUE coordinates, never scaled

  # ---- minimum drawn width --------------------------------------------------
  # A 2 kb gene in a 1.5 Mb window is a third of a pixel wide, so it is simply
  # absent from the figure. Growing only its exons does not help: the whole
  # model is that narrow. Instead the feature is ZOOMED about its own midpoint
  # until it reaches the floor, and every part of it -- exons, CDS, chevrons,
  # arrowhead -- is mapped through the same factor, so a small gene reads as a
  # small gene rather than as a smear.
  #
  # This is the one place the drawing stops being to scale, so it is opt-in and
  # anchored on the midpoint: the feature stays exactly where it belongs and
  # only its span is exaggerated.
  min_fw <- max(a$min_feature_width, a$min_feature_frac * window)
  fmid   <- (tstart + tend) / 2
  wf     <- widen(tstart, tend, min_fw)
  # how much each feature was stretched; 1 = untouched
  span   <- tend - tstart
  zoomf  <- ifelse(span > 0, (wf$xmax - wf$xmin) / span, 1)
  # map any coordinate inside feature k onto the zoomed model
  zoom   <- function(x, k) fmid[k] + (x - fmid[k]) * zoomf[k]
  fstart <- wf$xmin; fend <- wf$xmax

  # ---- which features get a label -------------------------------------------
  # Chosen BEFORE the vertical layout, because the answer decides how much of
  # the band the labels are entitled to. Selection reads the TRUE widths: a
  # gene does not become worth naming because it was drawn wider.
  lab_cfg <- a$labels
  cand <- data.table(idx = seq_len(n_feat), start = tstart, end = tend,
                     gene = nms, score = scr)[!is.na(gene) & nzchar(gene)]
  sel  <- select_labels(cand, lab_cfg, window)

  # ---- vertical placement inside this track's band --------------------------
  # `packed` splits the band into as many rows as the features need; `single`
  # and `jitter` keep one row. Feature height is always relative to the ROW, so
  # adding rows makes each glyph proportionally thinner instead of overflowing.
  # Packing reads the DRAWN extent, so two hairline genes blown up to the
  # minimum width are still kept off each other's row.
  if (a$layout == "packed") {
    rows   <- pack_rows(fstart, fend, a$row_gap_frac * window, a$max_rows)
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
  #
  # WHERE THE LABELS LIVE -- two independent choices:
  #
  #   anchor: band  one strip for the whole track, at one end of the band
  #           row   every label rides its OWN row, just clear of its own feature
  #   side:   top   the strip, or the rise, goes above the features
  #           bottom          ...below them
  #
  # `band` is the safe default: features and labels get disjoint zones, so no
  # label can ever land on a gene. Its cost is that the strip is at one END of
  # the band, so with packed rows the far row's labels are a whole band away
  # from what they name, and no amount of repulsion force can pull them closer
  # -- the strip is exactly where they were told to be. `row` trades that
  # guarantee for short leader lines: a label sits one row-gap from its own
  # feature, whichever row that is.
  anchor <- tolower(as.character(lab_cfg$anchor %||% "band"))
  side   <- tolower(as.character(lab_cfg$side   %||% "top"))
  if (!anchor %in% c("band", "row")) {
    stop("tracks.list.", a$name, ": 'labels.anchor: ", anchor,
         "' is not valid. Use one of: band, row", call. = FALSE)
  }
  if (!side %in% c("top", "bottom")) {
    stop("tracks.list.", a$name, ": 'labels.side: ", side,
         "' is not valid. Use one of: top, bottom", call. = FALSE)
  }

  # An empty strip is dead space, so the reserve is only taken when labels are
  # actually going to be drawn: turn `labels.enabled` off, or land on a locus
  # where every candidate was filtered out, and the genes get the whole band
  # instead of the fraction left over by a strip nothing will occupy. Under
  # `anchor: row` there is no strip at all -- the labels are among the features
  # -- so the reserve is not taken there either.
  reserve <- if (nrow(sel) && anchor == "band") {
               max(0, min(0.9, as.numeric(lab_cfg$reserve %||% 0)))
             } else 0
  feat_h  <- a$height * (1 - reserve)          # height left for the features
  # The features sit at the end of the band OPPOSITE the strip, so `feat_top`
  # is the top of the feature zone in both cases.
  feat_top <- if (side == "bottom") a$ymax else a$ymin + feat_h

  row_h  <- feat_h / n_rows                    # height of one row
  half   <- row_h / 2                          # half-height of a row
  h      <- half * a$feature_height            # half-height of a drawn feature
  utrh   <- if (isTRUE(a$utr_thinner)) h * a$utr_frac else h
  minw   <- a$min_exon_frac * window

  # row 1 sits at the TOP of the feature zone, like a genome browser
  row_centre <- function(r) feat_top - (r - 0.5) * row_h
  ybase <- row_centre(rows)
  if (a$layout == "jitter" && a$jitter > 0) {
    if (!is.null(a$seed)) set.seed(a$seed)
    ybase <- ybase + runif(n_feat, -half * a$jitter, half * a$jitter)
  }

  # Each label is nudged to the centre of a ZONE and repulsion may then move it
  # within that zone. `ylim` sets how much of the zone repulsion may use, and
  # `offset` slides the zone. What differs between the two anchors is how many
  # zones there are and where they sit.
  band_half <- a$height / 2
  gap    <- max(0, min(0.95, as.numeric(lab_cfg$gap %||% 0)))
  ylimf  <- as.numeric(lab_cfg$ylim %||% 1)
  offs   <- band_half * as.numeric(lab_cfg$offset %||% 0)
  rise_sign <- if (side == "bottom") -1 else 1

  if (anchor == "row") {
    # One zone per row, each just clear of that row's own features: past the
    # glyph's own half-height (feature_height) and `gap` row-halves further.
    # Neighbouring zones stay disjoint as long as ylim <= 1, so a label still
    # cannot drift into the row above or below.
    rise  <- half * (a$feature_height + gap)
    slack <- half * ylimf
    zone_of <- rows
    zone_ids <- sort(unique(rows))
    label_zones <- data.table(
      zone   = zone_ids,
      centre = row_centre(zone_ids) + rise_sign * rise + offs,
      slack  = slack)
  } else {
    # One strip for the whole track. `gap` raises its FLOOR (the edge nearest
    # the features), so no label can settle closer than that to the genes
    # however crowded the locus.
    if (reserve > 0) {
      strip_h <- a$height * reserve
      if (side == "bottom") {
        zone_lo <- a$ymin
        zone_hi <- (a$ymax - feat_h) - strip_h * gap
      } else {
        zone_lo <- feat_top + strip_h * gap
        zone_hi <- a$ymax
      }
    } else {
      zone_lo <- a$ycen - band_half
      zone_hi <- a$ycen + band_half
    }
    zone_of  <- rep(1L, n_feat)
    label_zones <- data.table(
      zone   = 1L,
      centre = (zone_lo + zone_hi) / 2 + offs,
      slack  = (zone_hi - zone_lo) / 2 * ylimf)
  }

  # `start`/`end` here are the DRAWN extent: everything below is geometry, and
  # the true coordinates have already done their job in the label selection.
  feat <- data.table(start = fstart, end = fend, strand = strd,
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
      # BED block coordinates are relative to the feature start, and are laid
      # out on the TRUE coordinates before being mapped onto the drawn model --
      # only then does an exon keep its position within a zoomed gene.
      xs <- zoom(tstart[k] + start(b) - 1, k)
      xe <- zoom(tstart[k] + end(b) - 1, k)
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
        c1 <- pmax(e$xmin, zoom(start(thick)[k], k))
        c2 <- pmin(e$xmax, zoom(end(thick)[k], k))
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
  # `sel` was chosen before the layout; all that is left is to put each label
  # where its feature ended up.
  lab <- lab_cfg
  labels <- empty
  if (nrow(sel)) {
    # `offset` is a fraction of the BAND half-height, not of the row: a label
    # is pushed away from the features as a group, and with packed rows a
    # row-relative offset would shrink with every extra row -- exactly when you
    # need it most. `jitter` scatters within the label's own zone.
    sel_zone <- zone_of[sel$idx]
    zsl <- label_zones[match(sel_zone, zone), slack]
    zc  <- label_zones[match(sel_zone, zone), centre]
    jt  <- zsl * as.numeric(lab$jitter %||% 0)
    if (!is.null(lab$seed)) set.seed(lab$seed)
    sel_y <- ybase[sel$idx]
    labels <- data.table(
      x = fmid[sel$idx],
      # The anchor stays ON the feature. ggrepel draws its leader line from
      # the label back to this point, so baking the offset into y here would
      # aim every line at empty space above the gene instead of at the gene.
      y = sel_y,
      # Nudge each label from its feature to the centre of its zone, so
      # repulsion starts where the label belongs instead of on top of the gene.
      nudge = (zc - sel_y) +
              (if (any(jt > 0)) runif(nrow(sel), -1, 1) * jt else 0),
      zone = sel_zone,
      label = sel$gene, track = a$name)
  }

  list(boxes = boxes, introns = introns, exons = exons, cds = cds,
       chevrons = chevrons, triangles = triangles, labels = labels,
       # computed once here so the renderer never re-derives the geometry
       label_zones = label_zones, score_lims = score_lims)
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
  window <- visible_region$end - visible_region$start

  # WHAT IS QUERIED IS NOT WHAT IS DRAWN. The visible interval may run past
  # either end of the chromosome -- `padding_fraction` adds context that does
  # not exist at a telomere, and a whole-chromosome locus is normally written
  # past the end -- and a BigWig query whose start is below 1 returns NOTHING,
  # with no error and no warning: the mark simply comes out blank. So the file
  # queries are clamped to the chromosome while `window` and the coordinates
  # keep the full visible range, which is what holds the tracks in step with
  # the heatmaps above.
  chr_size <- chr.list[chr == chr_id, size]
  q_start  <- max(1, visible_region$start)
  q_end    <- if (length(chr_size)) min(chr_size, visible_region$end) else visible_region$end
  if (q_end <= q_start) {
    stop("region '", locus$region, "' lies outside chr", chr_id,
         ", which is ", chr_size, " bp long", call. = FALSE)
  }
  region_gr <- GRanges(seqnames = chr_label, ranges = IRanges(q_start, q_end))

  # ---- marks ----------------------------------------------------------------
  mk_tracks <- tracks_of("mark")
  bw_signal_dt <- bin_marks(region_gr, locus)
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
      # ONE LAYER PER ZONE. ggrepel takes a single `ylim` for a whole layer, so
      # a track whose labels ride their own rows needs one layer per row --
      # otherwise every label would inherit one row's limits. The zones are
      # disjoint by construction, so splitting them costs no repulsion that
      # could have mattered: two labels that cannot reach the same height
      # cannot collide.
      zones <- an_parts[[i]]$label_zones
      for (zid in sort(unique(d$zone))) {
        dz <- d[zone == zid]
        zc <- zones[zone == zid]
        tracks_plot <- tracks_plot + geom_text_repel(
          data = dz, aes(x = x, y = y, label = label),
          # nudge_y moves the label while leaving the anchor on the feature, so
          # the leader line keeps pointing at the gene it names.
          nudge_y = dz$nudge,
          size = as.numeric(lab$size %||% 2),
          color = as.character(lab$color %||% "black"),
          fontface = as.character(lab$fontface %||% "bold"),
          # `force` is how hard labels shove EACH OTHER apart, `force_pull` how
          # hard each is drawn back to its own feature. Both only ever act
          # WITHIN the zone: `ylim` below is a hard boundary, so neither force
          # can pull a label out of the strip it was assigned. If the labels sit
          # too far from the genes, that is the zone's doing -- move it with
          # `anchor`, `side`, `reserve` and `gap`, not with the forces.
          force = as.numeric(lab$force %||% 15),
          force_pull = as.numeric(lab$force_pull %||% 1),
          max.overlaps = Inf,
          # direction: "x" keeps every label at its height and spreads them
          # sideways -- usually what you want above a gene track; "y" stacks
          # them vertically; "both" lets repulsion use either axis.
          direction = as.character(lab$direction %||% "both"),
          min.segment.length = as.numeric(lab$min_segment_length %||% 0),
          segment.size = as.numeric(lab$segment_size %||% 0.2),
          segment.color = as.character(lab$segment_color %||% lab$color %||% "black"),
          box.padding = as.numeric(lab$box_padding %||% 0.2), point.padding = 0.1,
          ylim = c(zc$centre - zc$slack, zc$centre + zc$slack), inherit.aes = FALSE)
      }
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

  # 6) colour bars for the annotation tracks drawn on a gradient. Added last:
  # each one turns the track panel into a patchwork, and a theme applied to a
  # patchwork lands on every plot inside it, insets included.
  for (i in seq_along(an_tracks)) {
    ins <- build_score_legend(an_tracks[[i]], an_parts[[i]]$score_lims, tracks_ylim)
    if (!is.null(ins)) tracks_plot <- tracks_plot + ins
  }

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
