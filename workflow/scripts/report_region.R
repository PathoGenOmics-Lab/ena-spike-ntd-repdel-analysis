log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(Rsamtools)
library(tidyverse)
library(ggpubr)
library(ape)

Sys.setlocale("LC_TIME", "English")
theme_set(theme_minimal())

NT.PALETTE <- c(
  "A" = "#118AB2",
  "C" = "#06D6A0",
  "G" = "#F78C6B",
  "T" = "#FFD166",
  "-" = "#073B4C",
  "+" = "#EF476F"
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

message("Reading reference FASTA file")
ref.seq <- read.dna(
  snakemake@input$reference,
  format = "fasta",
  as.character = TRUE)[snakemake@params$positions]

message("Building reference data")
reference <- data.frame(
  pos = seq_along(ref.seq) + snakemake@params$positions[1] - 1,
  ref = toupper(ref.seq)
)

message("Filtering and joining reference sequence")
plot.data <- bam.pileup %>%
    filter(pos %in% snakemake@params$positions) %>%
    left_join(reference, by = "pos")
    rename(count_raw = count) %>%
    group_by(pos) %>%
    mutate(count = ifelse(
        nucleotide == ref & "+" %in% nucleotide,
        count_raw - count_raw[nucleotide == "+"],
        count_raw
      )
    ) %>%
    ungroup()

message("Plotting")
p <- plot.data %>%
    ggplot(aes(pos, count, fill = nucleotide)) +
    geom_col(position = "stack") +
    geom_text(aes(label = ref, y = -30, color = ref), size = 3) +
    geom_text(aes(label = "Reference:", y = -30, x = snakemake@params$positions[1] - 5), size = 3) +
    scale_fill_manual(values = NT.PALETTE) +
    scale_color_manual(values = NT.PALETTE) +
    scale_x_continuous(breaks = snakemake@params$positions) +
      theme(
      axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5),
      panel.grid.minor = element_blank()
    )

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
