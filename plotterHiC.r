#!/usr/bin/env Rscript
# plotterHiC.r
# A minimalist R tool for visualizing differential Hi-C matrices in a 45° rotated square layout together with epigenomic tracks and gene annotations.
# Author: Oscar Amaury Aguilar Lomas
# Date: May 5, 2026

# Usage:
# Rscript plotterHiC.r configuration-file region-to-plot prefix_in_output

# Example
# Rscript plotterHiC.r plotter_config.sh chr1:78000000-78600000 Tst

# Bibliotecas necesarias -----------
library(dplyr)
library(ggplot2)
library(data.table)
library(parallel)
library(grid)
library(raster)
library(patchwork)
library(GenomicRanges)
library(yaml)
library(rtracklayer)
#
variables <- commandArgs(trailingOnly = TRUE)
config <- yaml::read_yaml(as.character(variables[1]))
#config <- yaml::read_yaml("config/plotter_config.yaml")

locus_coordinates <- strsplit(variables[2], ":")
(chr <- locus_coordinates[[1]][1])
chr <- gsub("chr", "", chr)
(start <- strsplit(locus_coordinates[[1]][2],"-")[[1]][1])
(end <- strsplit(locus_coordinates[[1]][2],"-")[[1]][2])

nomen <- variables[3]
# Si no está instalado, este es el comando para instalar straw :
# remotes::install_github("aidenlab/straw/R")
mm10 <- data.table(chr = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y"), size = c(195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271, 90702639, 61431566, 171031299, 91744698))
mm10_test <- data.table(chr = c("18", "19"), size = c(90702639, 61431566))
hg38 <- data.table(chr = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"), size = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415))
hg19 <- data.table(chr = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"), size = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566))
genome.list <- list(mm10 = mm10, hg38 = hg38, hg19 = hg19, mm10_test = mm10_test)

#------------------------------ Variables --------------------------------------
#General
(chr.list.var <- config$genome_reference)
# HiC
(hic.dir <- config$hic$directory)
(hic.file.ctrl <- config$hic$ctrl_file)
(hic.file <- config$hic$file)
(hic_names <- config$hic$names)
(matrix_choices <- config$hic$matrix_choices)
(resol <- config$hic$resolution)
(min_plot <- as.numeric(config$hic$min_plot))
(max_plot <- as.numeric(config$hic$max_plot))
# Heatmap legends
(show_legend <- as.logical(config$heatmap_legends$show))
(legend_position_x <- as.numeric(config$heatmap_legends$position_x))
(legend_position_y <- as.numeric(config$heatmap_legends$position_y))
(legend_direction <- config$heatmap_legends$direction)
(legend_title_size <- as.numeric(config$heatmap_legends$title_size))
(legend_text_size <- as.numeric(config$heatmap_legends$text_size))
(legend_key_height <- as.numeric(config$heatmap_legends$key_height))
(legend_key_width <- as.numeric(config$heatmap_legends$key_width))
# Marks
(marks_dir <- config$marks$directory)
(marks_names <- config$marks$names)
(marks <- config$marks$files)
(marks_colors <- config$marks$colors)
(track_offset <- as.numeric(config$marks$track_offset))
(track_height <- as.numeric(config$marks$track_height))
(marks_label_size <- config$marks$label_size)
(marks_bin_size <- as.numeric(config$marks$bin_size))
# Genes
(genes_dir <- config$genes$directories)
(genes_files <- config$genes$genes_files)
(genes_offset <- as.numeric(config$genes$offset))
(genes_labels <- config$genes$labels)
(genes_labels_offset <- as.numeric(config$genes$labels_offset))
(gene_height <- as.numeric(config$genes$height))
(gene_label_size <- as.numeric(config$genes$label_size))
# Plot
(plot_proportions <- as.numeric(config$plot$proportions))
(dir_plots <- config$plot$directory)
(plot_format <- config$plot$format)
(plot_height <- config$plot$height)
(plot_width <- config$plot$width)
# Shared groups
(shared_groups <- config$shared_yaxis_groups)

mat_parameters <- data.table(
   hic_names=hic_names,
   choices=matrix_choices,
   min_plot=min_plot,
   max_plot=max_plot,
   show_legend=show_legend,
   legend_position_x=legend_position_x,
   legend_position_y=legend_position_y,
   legend_direction=legend_direction,
   legend_title_size=legend_title_size,
   legend_text_size=legend_text_size,
   legend_key_height=legend_key_height,
   legend_key_width=legend_key_width)

chr.list <- genome.list[[chr.list.var]]
#------------------------------- Measure Time --------------------------------
start.time <- Sys.time()
options(scipen = 999)
dir <- getwd()
setwd(dir)
#------------------------------- FUNCTIONS -----------------------------------
# Generacion de polygonos -----------------
to_polygon <- function(x, y, counts) {
    # Four corners before rotation
    half <- (resol / 2)
    corners <- data.frame(
        x = c(x - half, x + half, x + half, x - half),
        y = c(y - half, y - half, y + half, y + half)
    )
    # Rotation of -45 degrees
    rot <- function(x, y) {
        u <- (x + y) / sqrt(2)
        v <- (y - x) / sqrt(2)
        return(data.frame(u = u, v = v))
    }
    rotated <- rot(corners$x, corners$y)
    rotated$counts <- counts
    rotated$id <- paste(x, y, sep = "_")
    rotated$corner <- 1:4
    rotated
}
# colores
scale_colors <- function(hicName, dt, values, minimo = NA, maximo = NA) {
    if (is.numeric(minimo) & is.numeric(maximo)) {
        if (values == "counts" | values == "counts_ctrl" | values == "log2" | values == "log2_ctrl") {
            scale <- scale_fill_gradientn(
                name = as.character(paste0(hicName,"\nInteractions \n(log2)")),
                colours = c("white", "orange", "darkred", "black"),
                limits = c(minimo, maximo),
                oob = scales::censor,
                na.value = "transparent"
            )
        } else if (values == "subtraction") {
            scale <- scale_fill_gradientn(
                name = paste0(as.character("subtraction\n",hicName)),
                limits = c(minimo, maximo),
                colours = c("blue", "white", "red"),
                values = scales::rescale(c(minimo, 0, maximo))
            )
        } else if (values == "ratio") {
            scale <- scale_fill_gradientn(
                name = paste0(as.character("Ratio",hicName)),
                limits = c(minimo, maximo),
                colours = c("blue", "white", "red"),
                values = scales::rescale(c(minimo, 1, maximo))
            )
        }
    } else {
        if (values == "counts" | values == "counts_ctrl" | values == "log2" | values == "log2_ctrl") {
            scale <- scale_fill_gradientn(
              as.character(paste0(hicName,"\nInteractions \n(log2)")),
                colours = c("white", "orange", "darkred", "black")
            )
        } else if (values == "subtraction") {
            scale <- scale_fill_gradientn(
              name = paste0(as.character("subtraction\n",hicName)), 
                limits = c(min(dt$counts, na.rm = TRUE), max(dt$counts, na.rm = TRUE)),
                colours = c("blue", "white", "red"),
                values = scales::rescale(c(min(dt$counts, na.rm = TRUE), 0, max(dt$counts, na.rm = TRUE)))
            )
        } else if (values == "ratio") {
            scale <- scale_fill_gradientn(
              name = paste0(as.character("Ratio",hicName)), 
                limits = c(min(dt$counts, na.rm = TRUE), max(dt$counts, na.rm = TRUE)),
                colours = c("blue", "white", "red"),
                values = scales::rescale(c(min(dt$counts, na.rm = TRUE), 1, max(dt$counts, na.rm = TRUE)))
            )
        }
    }
    return(scale)
}
# Funcion para invertir y graficar matriz horizontal
matriz_horizontal <- function(hicName,
                              plot_df,
                              region,
                              values, 
                              minimo = NA, 
                              maximo = NA,
                              show_legend = TRUE,
                              legend_position = c(0, 0.5),
                              legend_direction = "vertical",
                              legend_title_size = 8,
                              legend_text_size = 7,
                              legend_key_height = 20,
                              legend_key_width = 8) {
    plot_df <- plot_df[, .(x, y, counts = get(values))]
    plot_df$x <- as.numeric(as.character(plot_df$x))
    plot_df$y <- as.numeric(as.character(plot_df$y))
    plot_df$u <- (plot_df$x + plot_df$y) / sqrt(2)
    plot_df$v <- (plot_df$y - plot_df$x) / sqrt(2)
    # Toma el triangulo superior de la matriz
    plot_df_tri <- plot_df[plot_df$y >= plot_df$x, ]
    # Convierte a poligonos
    polygonos <- plot_df_tri[, to_polygon(x, y, counts), by = .(x, y)]
    # Coordenadas originales
    quantil_v <- quantile((unique(plot_df$v)), c(0, .7))
    quantil_u <- quantile((unique(plot_df$u)), c(.2, .8))
    
    # Gradiente de colores
    scale <- scale_colors(hicName,polygonos, values, minimo, maximo)

    # Plot
    apa_plot_h <- ggplot(polygonos, aes(u, v, group = id, fill = counts)) +
        geom_polygon() +
        coord_cartesian(clip = "off") +
        scale +
        ggtitle(region) +
        scale_y_continuous(limits = c(0, quantil_v[2]), expand = c(0, 0)) +
        scale_x_continuous(limits = c(quantil_u[1], quantil_u[2]),
                           position = "top",expand = c(0, 0),
                           labels = function(u_vals) {round(u_vals / sqrt(2) / 1e6, 1)}) +
        theme_minimal() +
        theme(
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),
            axis.title.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 10),
            axis.line = element_blank(),
            axis.ticks.x = element_line(color = "black", linewidth = 0.5),
            axis.ticks.length.x = unit(3, "pt"),
            axis.ticks.y = element_blank(),
            axis.ticks.length.y = element_blank(),
            panel.border = element_rect(linetype = "solid", fill = NA),
            plot.margin = margin(0,0,0,20),
            legend.position = if(show_legend) legend_position else "none",
            legend.direction = legend_direction,
            legend.title = element_text(size = legend_title_size,
                                        face = "bold",hjust = 0.5),
            legend.text = element_text(size = legend_text_size,
                                       face = "bold",hjust = 0.5),
            legend.key.height = unit(legend_key_height, "pt"),
            legend.key.width = unit(legend_key_width, "pt")
        )
    
    visible_coords <- list(
    start = quantil_u[1] / sqrt(2),
    end   = quantil_u[2] / sqrt(2))
    
    return(list(
    plot = apa_plot_h,
    visible = visible_coords))
}

u_to_genomic <- function(u){
    u / sqrt(2)
}
genomic_to_u <- function(x){
    x * sqrt(2)
}
#-----------------------------------------------------------------------------
chr_idx <- chr
bins.chr <- (round(chr.list[as.integer(chr_idx), 2] / resol)**2) / 2

locus <- c(78000000, 78600000)
print(paste0("Cargando matrices .hic de chr", chr))
mtz <- strawr::straw("KR", paste0(hic.dir, hic.file), chr, chr, "BP", resol, matrix = "observed")
mtz_ctrl <- strawr::straw("KR", paste0(hic.dir, hic.file.ctrl), chr, chr, "BP", resol, matrix = "observed")

avg_o <- sum(mtz$counts,na.rm=T) / bins.chr$size
avg_o_ctrl <- sum(mtz_ctrl$counts,na.rm=T) / bins.chr$size
avg_sum <- (avg_o / 2) + (avg_o_ctrl / 2)

setDT(mtz)
setDT(mtz_ctrl)

mtz[, norm := (counts / avg_o)]
mtz_ctrl[, norm := (counts / avg_o_ctrl)]
mtz[, log2 := log2(norm)]
mtz_ctrl[, log2 := log2(norm)]
mtz[, counts := NULL]
mtz_ctrl[, counts := NULL]
setnames(mtz, "norm", "counts")
setnames(mtz_ctrl, "norm", "counts")

setkey(mtz, x, y)
setkey(mtz_ctrl, x, y)
mtz_diff <- merge(mtz, mtz_ctrl, by = c("x", "y"), all.x = TRUE, suffixes = c("", "_ctrl"))
mtz_diff[, subtraction := (counts - counts_ctrl) * avg_sum]
mtz_diff[, ratio := counts / counts_ctrl]

anchor_down <- seq.int(start, end, resol)
anchor_up <- seq.int(start, end, resol)
cartesian_map <- CJ(x = anchor_down, y = anchor_up, sorted = FALSE)
setkey(cartesian_map, x, y)
locus_map <- mtz_diff[cartesian_map]
num_cols <- names(locus_map)[sapply(locus_map, is.numeric)]
locus_map[, (num_cols) := lapply(.SD, nafill, fill = 0), .SDcols = num_cols]

#-------------------------------------------------------------------------------
mat_ls <- list()
for(i in seq_along(matrix_choices)){
  mat_ls[[i]] <- 
    matriz_horizontal(
      mat_parameters$hic_names[i],
      locus_map,
      paste0("chr", chr),
      mat_parameters$choices[i],
      mat_parameters$min_plot[i],
       mat_parameters$max_plot[i],
       mat_parameters$show_legend[i],
       c(mat_parameters$legend_position_x[i], 
         mat_parameters$legend_position_y[i]),
       mat_parameters$legend_direction[i],
       mat_parameters$legend_title_size[i],
       mat_parameters$legend_text_size[i],
       mat_parameters$legend_key_height[i],
       mat_parameters$legend_key_width[i])
}

plot_ls <- list()
for(i in seq_along(matrix_choices)){
  if(i == 1){
    plot_ls[[i]] <- mat_ls[[i]]$plot +
      theme(
        axis.text.x.top = element_text(size = 6,
                                       margin = margin(b = -1)),
        axis.ticks.x.top = element_line(color = "black"),
        plot.margin = margin(0,0,0,20),
        plot.title = element_text(size = 6,
                                  margin = margin(b = -1))
      )
  } else {
  plot_ls[[i]] <- mat_ls[[i]]$plot +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        plot.margin = margin(0,0,0,0)
    )
  }
}

# visible genomic interval from HiC
visible_region <- mat_ls[[1]]$visible

#----------------------------- IMPORT BIGWIGS ----------------------------------
bw_files <- as.list(marks)
names(bw_files) <- marks_names
names(marks_colors) <- marks_names
names(track_offset) <- marks_names
names(track_height) <- marks_names
names(marks_bin_size) <- marks_names

region_gr <- GRanges(
    seqnames = paste0("chr", chr),
    ranges = IRanges(
        start = visible_region$start,
        end   = visible_region$end
    )
)
#----------------------------- BIN SIGNAL -----------------------------

bw_signal_ls <- list()
# bins <- tile(region_gr, n = nbins)[[1]]
for(mark in names(bw_files)){
  
  bins <- tile(region_gr, n = marks_bin_size[mark])[[1]]

    bw <- import(
        paste0(marks_dir, bw_files[[mark]]),
        which = region_gr
    )

    signal_vec <- numeric(length(bins))

    for(i in seq_along(bins)){

        overlaps <- subsetByOverlaps(
            bw,
            bins[i]
        )

        signal_vec[i] <- ifelse(
            length(overlaps) > 0,
            mean(overlaps$score),
            0
        )
    }

    bw_dt <- data.table(
        start = start(bins),
        end   = end(bins),
        signal = signal_vec,
        mark = mark
    )

    bw_signal_ls[[mark]] <- bw_dt
}

bw_signal_dt <- rbindlist(bw_signal_ls)

#----------------------------- TRACK CONFIG ------------------------------------
track_config <- data.table(
    mark = names(bw_files),
    y_offset = as.vector(track_offset[names(bw_files)]),
    track_height = track_height[names(bw_files)],
    color = marks_colors[names(bw_files)],
    label_offset = 0.4
)
#----------------------------- LOG SIGNAL --------------------------------------
bw_signal_dt[, log_signal := log10(1 + signal)]
#----------------------------------Share scale----------------------------------
for(group in shared_groups){
  shared_max <- bw_signal_dt[
    mark %in% group,
    max(log_signal, na.rm = TRUE)
  ]
  bw_signal_dt[
    mark %in% group,
    display_max := shared_max
  ]
}
# merge configuration
setDT(bw_signal_dt)
setDT(track_config)
bw_signal_dt <- track_config[bw_signal_dt, on = "mark"]
# independent scaling per track
bw_signal_dt[, ymax := y_offset + (log_signal * track_height)]
bw_signal_dt[, ymin := y_offset]
#----------------------------- GENE TRACK --------------------------------------

bed_file <- "/Users/amaury/genes/Ctcf_WT_KO.DEGs.filtered_by_Tie2_ATAC.Accesible_Inaccesible.genebody.bed"

genes_gr <- import(bed_file)

genes_gr <- genes_gr[
    seqnames(genes_gr) == paste0("chr", chr)
]

genes_gr <- subsetByOverlaps(
    genes_gr,
    region_gr
)

genes_dt <- data.table(
    start = start(genes_gr),
    end   = end(genes_gr),
    strand = as.character(strand(genes_gr)),
    type = genes_gr$type,
    gene = genes_gr@elementMetadata$name,
    color = genes_gr@elementMetadata$itemRgb    
)

gene_track_y <- genes_offset
genes_dt[, y := gene_track_y]

# gene labels
gene_labels <- unique(
    genes_dt[, .(
        gene,
        start,
        end
    )]
)
gene_labels[, middle := (start + end)/2]

gene_arrows <- unique(
    genes_dt[, .(
        gene,
        strand,
        start,
        end,
        color
    )]
)

gene_arrows[, arrow_x := ifelse(
    strand == "+",
    end,
    start
)]

gene_arrows[, arrow_dir := ifelse(
    strand == "+",
    1,
    -1
)]

gene_arrows[, y := gene_track_y]
#-------------------------
gene_box_height <- gene_height
arrow_height <- gene_box_height*3
arrow_ymin <- -gene_box_height
arrow_ymax <-  gene_box_height

triangle_dt <- gene_arrows[, {
    if(strand[1] == "+"){
        data.table(
            x = c(
                end[1],
                end[1] + 6000,
                end[1]
            ),
            y = c(
                y[1] + arrow_ymin,
                y[1] + (arrow_ymin + arrow_ymax)/2,
                y[1] + arrow_ymax
            ), color = color[1]
        )
    } else {
        data.table(
            x = c(
                start[1],
                start[1] - 6000,
                start[1]
            ),
            y = c(
                y[1] + arrow_ymin,
                y[1] + (arrow_ymin + arrow_ymax)/2,
                y[1] + arrow_ymax
            ),color = color[1]
        )
    }
}, by = gene]
#----------------------------- BIGWIG BROWSER PLOT -----------------------------

# track_labels <- bw_signal_dt[
#   ,.(max=round(max(log_signal),2),
#      y_offset=unique(y_offset)),by = mark]

track_labels <- bw_signal_dt[
  ,.(max = round(unique(display_max),2),
     y_offset = unique(y_offset)),
  by = mark]

track_distance <- visible_region$end-visible_region$start
strat_track <- visible_region$start
genes_labels_dt <- data.table(mark = genes_labels, y_offset = genes_labels_offset)
track_labels_dt <- rbind(track_labels[,.(mark,y_offset)],genes_labels_dt)
tracks_plot <- ggplot() +
    # signal tracks
    geom_ribbon(
    data = bw_signal_dt,
    aes(
        x = start,
        ymin = ymin,
        ymax = ymax,
        fill = color,
        group = mark
    ),
    alpha = 0.9,
    color = "black",
    linewidth =0.2)+
    scale_fill_identity() +
    geom_rect(
        data = genes_dt,
        aes(
            xmin = start,
            xmax = end,
            ymin = y - gene_box_height,
            ymax = y + gene_box_height,
            fill = color
        ), color = "black") +
    geom_polygon(
    data = triangle_dt,
    aes(
        x = x,
        y = y,
        group = gene,
        fill = color
    ),   
    color = "black",
    inherit.aes = FALSE)  +
    geom_text(data = gene_labels,
        aes(x = middle, y = genes_labels_offset, label = gene),
        size = gene_label_size,fontface = "bold") +
    scale_x_continuous(limits = c(visible_region$start, visible_region$end),
        expand = c(0,0),labels = function(x){round(x / 1e6, 2)},
        position = "top") +
  scale_y_continuous(limits = c(genes_labels_offset-0.5,max(track_offset)+1),
                     expand = c(0,0),
                     breaks = track_labels_dt$y_offset,
                     labels = track_labels_dt$mark)+
    theme_minimal() +
    theme(
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = marks_label_size,
                                   angle = 0,
                                   hjust = 1,
                                   face = "bold"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
        plot.margin = margin(0,0,0,1),
        panel.border = element_rect(
            fill = NA,
            color = "black"            
        )
    ) +
  geom_text(
    data = track_labels,
    aes(
      x = strat_track+(track_distance*0.005),
      y = y_offset + 0.85,
      label = max),
    inherit.aes = FALSE,
    size = 2,
    hjust = 0
  )+
  geom_text(
    data = track_labels,
    aes(
      x = strat_track+(track_distance*0.005),
      y = y_offset + 0.15,
      label = 0),
    inherit.aes = FALSE,
    size = 2,
    hjust = 0
  )
#----------------------------- COMBINE -----------------------------------------
final_plot <-
  plot_ls[[1]] / plot_ls[[2]] / plot_ls[[3]] /
    tracks_plot +
    plot_layout(
        heights = plot_proportions
    )

dir.create(dir_plots)

formato <- "png" # "pdf" o "png"
  if(formato=="png"){
    path.plot <- paste0(dir_plots,"/",nomen,".",
                        chr,"_",start,"_",end,".",
                        resol,".", "_bins.hicplot.png")
    device.plot <- png(path.plot,width=10,
                       height=15,units="cm",res=500)
  } else if (formato=="pdf"){
    path.plot <- paste0(dir_plots,"/",nomen,".",
                        chr,"_",start,"_",end,".",
                        resol,".", "_bins.hicplot.pdf")
    device.plot <- pdf(path.plot, width = 5 / 2.54, 
                       height = 10 / 2.54) 
  } else{
    stop("Image format not available. Use instead 'pdf' or 'png'")
  } 
  path.plot
  device.plot
  print(final_plot)
  dev.off()


