# plotterHiC

A minimalist R tool for visualizing Hi-C matrices in a 45° rotated layout, stacked over a genome-browser panel carrying bigwig signal, TAD/loop arcs and any BED annotation.

Everything is driven by one YAML file. Nothing about a figure lives in the script.

---

## Features

- One rotated heatmap per panel: raw counts, log2, subtraction or ratio
- A **named registry** of `.hic` files, so any file can be compared against any other
- Per-file resolution and normalization (files never have to agree)
- Colour limits as fixed numbers, as **quantiles** (`q0.01`), or automatic
- **Shared colour scales** across panels, so heatmaps are directly comparable
- One ordered track list: bigwig ribbons, TAD/loop arcs, gene models, any BED
- Full gene models (exons, introns, thinner UTRs, chevrons) or plain boxes
- Genome-browser **row packing** so overlapping features never collide
- Gene labels with isoform de-duplication, caps, and a reserved strip that keeps them off the genes
- `enabled: true/false` on every region, panel and track
- Automatic vertical layout: declaration order is figure order

---

## Example

<p align="center">
  <img src="figures/example_tracks.png" width="800">
</p>

Hi-C panels stacked over the track panel: bigwig signal, gene annotations and their labels, all aligned to the same genomic interval.

---

## Installation

```bash
git clone https://github.com/YOUR_USERNAME/plotterHiC.git
cd plotterHiC
```

### Dependencies

```r
install.packages(c("data.table", "ggplot2", "patchwork", "ggrepel", "yaml"))

install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "rtracklayer"))

# reads .hic files from Juicer
remotes::install_github("aidenlab/straw/R")
```

---

## Quick start

Two arguments, and you get a figure:

```bash
Rscript plotterHiC.r config/plotterHiC_mESC_ysEC.yaml chr2:112225000-114225000 Actc1
```

- `config/...yaml` — the configuration file
- `chr2:112225000-114225000` — the region (`chr2:112225000:114225000` also works)
- `Actc1` — a name for the output file (optional)

Or drop the region and let the config decide, drawing one figure per entry in `global.loci`:

```bash
Rscript plotterHiC.r config/plotterHiC_mESC_ysEC.yaml
```

To get going, copy an existing config and change three things: the paths under `hic.files`, the regions under `global.loci`, and the output `dir`.

---

## Input files

| Input | Used as |
|---|---|
| `.hic` (Juicer) | the heatmap panels |
| `.bw` (bigWig) | `type: mark` signal tracks |
| `.bedpe` | `type: structure` — loop arcs, TAD domains |
| `.bed` | `type: annotation` — gene models or plain boxes |

Genomes: `mm10`, `hg38`, `hg19`.

---

## ⚠ Write configs in ASCII only

`Rscript` launched from a terminal with `LANG` unset runs under the **C locale**. There, a single non-ASCII byte — an accented word in a *comment* is enough — makes the YAML reader truncate the file at that point and return only the keys above it, **with no error**.

plotterHiC reads config bytes directly to be immune to this, but keep configs ASCII anyway so they stay portable. Check with:

```bash
LC_ALL=C grep -n '[^ -~]' config/your_config.yaml
```

---

## Configuration

### How the file is shaped

Every panel and every track is a **braced block**, so where one ends and the next begins is unambiguous. Note the commas between entries, and that values containing `#` must be quoted (`"#1a1a1a"`).

```yaml
CTCF_mESC: {
  enabled: true,
  type: mark,
  file: /full/path/to/signal.bw,
  color: blue2
}
```

Three rules hold throughout:

- **Order is declaration order.** There are no y coordinates anywhere. Panels are drawn top to bottom as written; tracks likewise. To move something, move its block.
- **Every path is a full, independent path.** Nothing is assembled from a shared directory plus a name.
- **Whatever a block does not mention** falls back to the script's built-in default.

### The four sections

```yaml
global:  { genome, resolution, loci, output }
hic:     { norm, files, panel_defaults, legend_defaults, panels }
tracks:  { panel_height, spacing, label_size, list }
plot:    { x_axis, crop }
```

### `global.loci` — one figure per region

A locus is either a bare region or a block that overrides settings for that region alone:

```yaml
loci: {
  Actc1: chr2:112225000-114225000,          # inherits everything

  TADs_4Mb: {
    enabled: true,
    region: chr2:111000000-115000000,
    resolution: 25000, plot_width: 10, plot_height: 10,
    log2:        { min: q0.01, max: 13 },   # per matrix type
    subtraction: { min: q0.01, max: q1 }
  }
}
```

Locus keys: `enabled`, `region`, `resolution`, `padding_fraction`, `plot_width`, `plot_height`, `formato`, `dpi`, `dir`.

It may also carry panel overrides applied to every heatmap of that figure — `min`, `max`, `matrix`, `normalize`, `colors`, `midpoint`, `oob`, `na_color`, `height`, `title`, `legend` — either flat or inside a per-matrix block as above. A 4 Mb TAD view and a 30 Mb compartment view rarely share a colour range, and a log2 map and a rescaled subtraction never do.

`enabled: false` parks a region so a rerun regenerates only what changed. A region given on the command line is always drawn.

### `hic` — a registry, not a control/treatment pair

```yaml
hic: {
  norm: KR,

  files: {
    ysEC: /full/path/Ctcf_WT.hic,
    mESC: { file: /full/path/mESC.hic, resolution: 5000, norm: SCALE }
  },

  panel_defaults: {
    matrix: log2,
    log2:        { normalize: true },
    subtraction: { normalize: false }
  },

  panels: [
    { enabled: true,  name: mESC,          hic: mESC },
    { enabled: true,  name: ysEC,          hic: ysEC },
    { enabled: true,  name: "ysEC - mESC", hic: ysEC, vs: mESC, matrix: subtraction }
  ]
}
```

Panels reference files **by name**, so N files support any pairwise comparison. A comparison always reads as `hic <op> vs` — swap the two keys to flip it.

| `matrix` | Meaning | Needs |
|---|---|---|
| `counts` | contacts ÷ chromosome-wide mean | `hic` |
| `log2` | log2 of the above | `hic` |
| `subtraction` | `hic - vs`, rescaled to contact units | `vs` |
| `ratio` | `hic / vs` | `vs` |

Counts are normalized per file against that file's own chromosome-wide mean contact density, which is what makes libraries of different depth comparable.

Each file may pin its own `resolution` and `norm`. The one constraint: a `subtraction` or `ratio` panel merges two matrices bin by bin, so **both sides must share a resolution** — the script says so explicitly rather than drawing a misaligned map.

### Colour limits

`min` and `max` are read independently, so mixing forms is fine (`min: q0.01, max: 12`).

| Form | Meaning |
|---|---|
| `NA` | the data's own min / max in the visible window |
| a number | that exact value |
| `qN` | the Nth quantile, `N` between 0 and 1 |

Quantiles are the robust option: `q0.01` / `q0.99` clip the extreme 1% at each end, so one saturated bin cannot flatten the map, and no retuning is needed per locus.

**`normalize`** makes panels comparable. Panels in the same group take their limits from the **pooled** values of the whole group, so `NA` means the group's min/max and `qN` the group's quantile. `normalize: true` joins a default group; a string names one, for several independent groups in a figure.

### Where a setting can be written

Panel settings resolve in six layers, each narrower than the last:

```
hic.panel_defaults          (flat)
hic.panel_defaults[matrix]  per matrix type, figure-wide
hic.presets[preset]         optional named looks
the locus                   (flat)
the locus[matrix]           per matrix type, for this region
the panel                   always the last word
```

The panel stays last on purpose: a locus-wide `matrix: log2` must not silently retype a subtraction panel that pinned its own. The practical rule is **identity in the panel, scales in the locus** — a panel that pins `min`/`max` stops responding to its region.

`hic`, `vs`, `name`, `enabled` and `preset` cannot be set at locus level; they identify *which* panel a block is.

### `tracks` — one ordered list

```yaml
tracks: {
  panel_height: 2,     # the whole track panel vs one heatmap
  spacing: 0.2,        # gap between every consecutive pair
  label_size: 6,

  list: {
    Loops: { enabled: true, type: structure, file: /path/x.bedpe, style: arc,
             height: 1.2, max_n: 25, rank_by: length, rank_desc: true },

    CTCF:  { enabled: true, type: mark, file: /path/x.bw, color: blue2,
             height: 1, signal_scale: auto, space_above: 0.3 },

    RefSeq: { enabled: true, type: annotation, file: /path/genes.bed,
              style: model, layout: packed, max_rows: 4,
              direction: chevrons, feature_height: 0.6,
              labels: { enabled: true, dedup: true, max_n: 12,
                        reserve: 0.45, gap: 0.4 } }
  }
}
```

Each track owns a horizontal **band**. Three knobs, and nothing else, control the vertical rhythm:

| Knob | Scope | Effect |
|---|---|---|
| `height` | per track | how tall its band is; a mark ribbon fills its band |
| `spacing` | global | gap between every consecutive pair |
| `space_above` | per track | extra gap above this track only; negative pulls it up |

#### `type: mark` — bigwig signal

`signal_scale: auto` rescales the track to fill its band at every locus, which is what lets bigwigs from different assays sit side by side; a number fixes the multiplier instead. Also: `n_bins`, `crop`, `color`, `alpha`, `line_color`, `line_width`, `show_scale`.

Marks named together under `shared_yaxis_groups` share one y scale, so their heights are comparable.

#### `type: structure` — TADs and loops from a bedpe

`style: arc` draws a semicircle between anchor midpoints (loops); `style: domain` draws a bar spanning the domain (TADs). Filter with `max_n`, `rank_by` (`length` | `score` | `random`), `rank_desc`, `min_span`, `max_span`. Set `add_chr: true` when the bedpe stores `2` but the genome needs `chr2`. Non-standard column layouts are handled with `columns` and `select`.

#### `type: annotation` — any BED

`style: model` draws the full gene model from a BED12 — intron line, exon boxes, thinner UTRs, direction chevrons. `style: box` draws one rectangle per feature and works with any BED: peaks, motifs, enhancers.

`layout` decides how overlapping features are spread:

| `layout` | Behaviour |
|---|---|
| `single` | all on one line; overlaps draw on top of each other |
| `packed` | real rows — no two overlapping features share a row. Deterministic |
| `jitter` | random vertical scatter; shows density, but separation is luck |

With `packed`, feature height is relative to the **row**, so adding rows makes each glyph thinner rather than overflowing the band. `max_rows` caps them.

#### Gene labels

Labels are two separate problems, and the settings split the same way: **which** features get named, and **where** the names go.

```yaml
labels: {
  enabled: true,
  size: 1.5, color: black, fontface: italic,

  # WHICH features get a label
  dedup: true,            # one label per name — collapses isoforms
  max_n: 6,               # cap; 0 = no cap
  rank_by: length,        # length | score | position | random
  rank_desc: true,        # true = longest / highest first
  min_width_frac: 0.01,   # skip features under 1% of the visible window
  min_width: 0,           # ...or under this many bp
  include: [Nop10, Aven], # always label these, whatever the filters say
  exclude: [],            # never label these

  # WHERE the labels go
  reserve: 0.45,          # keep the top 45% of the band clear for labels
  gap: 0.4,               # floor: minimum clearance from the genes
  ylim: 0.9,              # how much of the remaining strip repel may use
  offset: 0,              # slide the whole label zone (band-half fraction)
  jitter: 0,              # extra random scatter within the zone
  repel: true,
  direction: both,        # both | x | y
  segment_size: 0.2,      # leader line to the feature; 0 removes it
  box_padding: 0.2
}
```

**Which features.** `dedup` runs before `include`, so a gene named in `include` is still labelled once rather than once per transcript. `include` bypasses the size filters and the `max_n` cap; everything else competes for the remaining slots by `rank_by`.

**Where they go.** This is what keeps labels off the genes:

| Setting | Effect |
|---|---|
| `reserve` | Splits the band: features pack into the bottom `1 - reserve`, labels live in the top `reserve`. The two zones are disjoint, so **no label can land on a gene**, however crowded the locus. `0` disables the split and labels share the band with the features. |
| `gap` | Raises the **floor** of the label zone, as a fraction of the reserved strip. This is the knob for "some labels end up stuck to the genes": repulsion is free to place a label anywhere it is allowed, so without a floor the crowded ones settle right on the features. |
| `ylim` | Width of the allowed interval. Smaller = more uniform label height, but less room to stagger and more crowding sideways. |
| `offset` | Slides the whole zone, in fractions of the band half-height. Useful with `reserve: 0` or to push labels out of the band entirely. |
| `direction` | `x` keeps every label at one height and spreads them sideways; `y` stacks them; `both` lets repulsion use either axis. |

Rule of thumb: for even, tidy heights raise `gap` and lower `ylim`; to let crowded labels stagger so they all fit, lower `gap` and raise `ylim`.

The leader line always points at the feature itself — the label is moved with a nudge while its anchor stays on the gene, so a line never aims at empty space.

### `plot`

```yaml
plot: {
  x_axis: { divide: 1000000, text_size: 6, title_size: 6 },
  crop:   { x: [0.2, 0.8], y: 0.7 }
}
```

`crop.x` trims the left and right corners of the rotated square, whose bins pair with positions outside the window; `crop.y` sets how tall a slice is kept. The crop is derived from the window itself, not from any panel's bins, which is what keeps panels of different resolutions exactly aligned.

---

## Turning things on and off

`enabled: false` works on **regions**, **panels** and **tracks**. In every case the item is skipped entirely — never validated, its file never opened — and whatever follows closes the gap. Keep everything you ever use in one config and flip a few booleans per figure.

---

## Output

```text
<dir>/<name>.<chr>_<start>_<end>.<resolution>._bins.hicplot.png
```

PNG or PDF, per figure, sized by `plot_width` / `plot_height` in cm.

---

## Repository structure

```text
plotterHiC/
├── plotterHiC.r
├── config/
│   ├── plotter_config.yaml            # Ctcf_WT vs Ctcf_KO
│   └── plotterHiC_mESC_ysEC.yaml      # mESC vs ysEC
└── README.md
```

---

## Author

Oscar Amaury Aguilar Lomas

## License

MIT License
