log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(Rsamtools)
library(tidyverse)
library(ggpubr)

Sys.setlocale("LC_TIME", "English")
theme_set(theme_minimal())

NT.PALETTE <- c(
  "A" = "#8b008b",
  "C" = "#ff7f50",
  "G" = "#6495ed",
  "T" = "#458b00",
  "-" = "#3d3d3d"
)

message("Setting params")
params <- PileupParam(
  max_depth = snakemake@params[["max_depth"]],
  min_base_quality = snakemake@params[["min_base_quality"]],
  min_mapq = snakemake@params[["min_mapq"]],
  min_nucleotide_depth = snakemake@params[["min_nucleotide_depth"]],
  min_minor_allele_depth = snakemake@params[["min_minor_allele_depth"]],
  distinguish_strands = snakemake@params[["distinguish_strands"]],
  distinguish_nucleotides = snakemake@params[["distinguish_nucleotides"]],
  ignore_query_Ns = snakemake@params[["ignore_query_Ns"]],
  include_deletions = snakemake@params[["include_deletions"]],
  include_insertions = snakemake@params[["include_insertions"]]
)

message("Reading BAM file")
bam <- BamFile(snakemake@input$bam)
bam.pileup <- pileup(bam, pileupParam = params)

message("Filtering")
plot.data <- bam.pileup %>%
    filter(pos %in% snakemake@param$positions)

message("Plotting")
p <- plot.data %>%
    ggplot(aes(pos, count, fill = nucleotide)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = NT.PALETTE)

message("Writing plot")
ggsave(
    snakemake@output$plot,
    plot = p,
    height = snakemake@params$plot_height_in,
    width = snakemake@params$plot_width_per_position_in *
        length(snakemake@params$positions),
    bg = "white"
)

message("Writing plot data")
write_csv(plot.data, snakemake@output$plot_data)
