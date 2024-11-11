log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)
library(ggpubr)

Sys.setlocale("LC_TIME", "English")
theme_set(theme_minimal())

message("Processing input tables")
plot.data <- lapply(
        snakemake@input$tables,
        function(path) {
            read_delim(path) %>%
                mutate(path = path) %>%
                filter(
                    POS %in% snakemake@params$positions,
                    PASS == snakemake@params$filter_pass
                )
        }
    ) %>%
    bind_rows() %>%
    distinct() %>%
    # Calculate the average ALT_FREQ across all samples
    group_by(REGION, POS, REF, ALT) %>%
    summarize(AVG_ALT_FREQ = mean(ALT_FREQ, na.rm = TRUE))

message("Plotting")
p <- plot.data %>%
    # Plot
    ggplot() +
    aes(x = as.factor(POS), y = AVG_ALT_FREQ, fill = ALT) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap("REGION")

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
