library(tidyverse)
library(ggpubr)

Sys.setlocale("LC_TIME", "English")
theme_set(theme_minimal())

format.date <- function(df = ., date.col, new.col) {
    # Dates with missing month get imputed as day 1 of that month
    date.col <- enquo(date.col)
    new.col <- enquo(new.col)
    df %>%
        mutate(
            (!!new.col) := case_when(
                str_count(!!date.col, "-") == 2 ~
                    as.Date(!!date.col, format = "%Y-%m-%d"),
                str_count(!!date.col, "-") == 1 ~
                    as.Date(paste0(!!date.col, "-01"), format = "%Y-%m-%d"),
                TRUE ~ NA
            )
        )
}

# Read search results
message("Reading data")
search <- read_tsv(
    snakemake@input$table,
    col_types = cols(
        .default = "c"
    ),
    n_max = 10000
)
message("Read ", nrow(search), " records")

# Impute missing collection dates
# DATE:
# collection_date, collection_date_end, collection_date_start
# first_public, last_updated
message("Imputing missing collection dates")
search <- search %>%
    format.date(collection_date, `Collection date`) %>%
    format.date(first_created, `First created date`) %>%
    format.date(collection_date_end, `Collection date end`) %>%
    format.date(collection_date_start, `Collection date start`) %>%
    format.date(last_updated, `Last updated`)

# Country timeline
# LOCATION:
# country
message("Plotting samples per month, by country")
date.country.data <- search %>%
    separate(country, into = c("Country", "Region"), sep = " *: *") %>%
    count(Country, Region, `Collection date`, `First created date`, `Collection date start`, `Collection date end`, `Last updated`) %>%
    pivot_longer(
        c(`Collection date`, `First created date`, `Collection date start`, `Collection date end`, `Last updated`),
        names_to = "DateType",
        values_to = "Date"
    ) %>%
    mutate(
        Location = paste(Country, Region, sep = "/"),
        Date = lubridate::floor_date(Date, unit = "month")
    )

date.country.p <- ggarrange(
    # Timeline
    date.country.data %>%
    ggplot() +
        aes(x = Date, y = n, fill = Country) +
        geom_bar(stat = "identity") +
        scale_x_date(date_breaks = "1 month", date_minor_breaks = "1 month") +
        facet_grid("DateType") +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5)
        ) +
        scale_fill_viridis_d() +
        guides(fill = guide_legend(ncol = 10)),
    # Missing dates
    date.country.data %>%
    mutate(`Has date` = !is.na(Date)) %>%
    ggplot() +
        aes(x = 1, fill = `Has date`) +
        geom_bar(stat = "count") +
        geom_text(
            aes(label = after_stat(count)),
            stat = "count",
            position = position_stack(vjust = 0.5)
        ) +
        facet_grid("DateType") +
        theme(
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        scale_fill_viridis_d(begin = 0.6, end = 0.9),
    # Missing country
    date.country.data %>%
    mutate(`Has country` = !is.na(Country)) %>%
    ggplot() +
        aes(x = 1, fill = `Has country`) +
        geom_bar(stat = "count") +
        geom_text(
            aes(label = after_stat(count)),
            stat = "count",
            position = position_stack(vjust = 0.5)
        ) +
        facet_grid("DateType") +
        theme(
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        scale_fill_viridis_d(begin = 0.6, end = 0.9),
    # Settings
    nrow = 1,
    align = "h",
    widths = c(8, 1, 1),
    legend = "top"
)

# Timeline by sequencing tech
# SEQUENCING:
# instrument_model, instrument_platform
# library_construction_protocol, library_gen_protocol,
# library_layout, library_name
# library_min_fragment_size, library_max_fragment_size
# library_pcr_isolation_protocol
# library_selection, library_source, library_strategy
# isolation_source
message("Plotting samples per month, by sequencing tech")
date.seqtech.data <- search %>%
    mutate(
        Tech = paste(instrument_platform, instrument_model, library_layout, library_strategy, sep = "-")
    ) %>%
    count(
        Tech,
        `Collection date`, `First created date`, `Collection date start`, `Collection date end`, `Last updated`
    ) %>%
    pivot_longer(
        c(`Collection date`, `First created date`, `Collection date start`, `Collection date end`, `Last updated`),
        names_to = "DateType",
        values_to = "Date"
    ) %>%
    mutate(
        Date = lubridate::floor_date(Date, unit = "month")
    )

date.seqtech.p <- ggarrange(
    # Timeline
    date.seqtech.data %>%
    ggplot() +
        aes(x = Date, y = n, fill = Tech) +
        geom_bar(stat = "identity") +
        scale_x_date(date_breaks = "1 month", date_minor_breaks = "1 month") +
        facet_grid("DateType") +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5)
        ) +
        scale_fill_viridis_d() +
        guides(fill = guide_legend(ncol = 3)),
    # Missing dates
    date.seqtech.data %>%
    mutate(`Has date` = !is.na(Date)) %>%
    ggplot() +
        aes(x = 1, fill = `Has date`) +
        geom_bar(stat = "count") +
        geom_text(
            aes(label = after_stat(count)),
            stat = "count",
            position = position_stack(vjust = 0.5)
        ) +
        facet_grid("DateType") +
        theme(
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        scale_fill_viridis_d(begin = 0.6, end = 0.9),
    # Missing sequencing tech
    date.seqtech.data %>%
    mutate(`Has tech` = !is.na(Tech)) %>%
    ggplot() +
        aes(x = 1, fill = `Has tech`) +
        geom_bar(stat = "count") +
        geom_text(
            aes(label = after_stat(count)),
            stat = "count",
            position = position_stack(vjust = 0.5)
        ) +
        facet_grid("DateType") +
        theme(
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        scale_fill_viridis_d(begin = 0.6, end = 0.9),
    # Settings
    nrow = 1,
    align = "h",
    widths = c(8, 1, 1),
    legend = "top"
)

# Timeline by sequencing results
# SEQUENCING:
# base_count, read_count
message("Plotting samples per month, by binned sequencing results")
date.seqres.data <- search %>%
    mutate(
        `Base count` = cut_interval(as.numeric(base_count), snakemake@params$count_bins),
        `Read count` = cut_interval(as.numeric(read_count), snakemake@params$count_bins)
    ) %>%
    count(
        `Base count`, `Read count`,
        `Collection date`, `First created date`, `Collection date start`, `Collection date end`, `Last updated`
    ) %>%
    pivot_longer(
        c(`Collection date`, `First created date`, `Collection date start`, `Collection date end`, `Last updated`),
        names_to = "DateType",
        values_to = "Date"
    ) %>%
    mutate(
        Date = lubridate::floor_date(Date, unit = "month")
    )

date.seqres.p <- ggarrange(
    # Timeline
    ggarrange(
        date.seqres.data %>%
        ggplot() +
            aes(x = Date, y = n, fill = `Base count`) +
            geom_bar(stat = "identity") +
            scale_x_date(date_breaks = "1 month", date_minor_breaks = "1 month") +
            facet_grid("DateType") +
            theme(
                axis.text.x = element_text(angle = 90, vjust = 0.5)
            ) +
            scale_fill_brewer(palette = "Set1") +
            guides(fill = guide_legend(ncol = 3)),
        date.seqres.data %>%
        ggplot() +
            aes(x = Date, y = n, fill = `Read count`) +
            geom_bar(stat = "identity") +
            scale_x_date(date_breaks = "1 month", date_minor_breaks = "1 month") +
            facet_grid("DateType") +
            theme(
                axis.text.x = element_text(angle = 90, vjust = 0.5)
            ) +
            scale_fill_brewer(palette = "Set1") +
            guides(fill = guide_legend(ncol = 3)),
        nrow = 1,
        align = "h",
        legend = "top"
    ),
    # Missing dates
    date.seqres.data %>%
    mutate(`Has date` = !is.na(Date)) %>%
    ggplot() +
        aes(x = 1, fill = `Has date`) +
        geom_bar(stat = "count") +
        geom_text(
            aes(label = after_stat(count)),
            stat = "count",
            position = position_stack(vjust = 0.5)
        ) +
        facet_grid("DateType") +
        theme(
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        scale_fill_viridis_d(begin = 0.6, end = 0.9),
    # Missing base counts
    date.seqres.data %>%
    mutate(`Has base n` = !is.na(`Base count`)) %>%
    ggplot() +
        aes(x = 1, fill = `Has base n`) +
        geom_bar(stat = "count") +
        geom_text(
            aes(label = after_stat(count)),
            stat = "count",
            position = position_stack(vjust = 0.5)
        ) +
        facet_grid("DateType") +
        theme(
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        scale_fill_viridis_d(begin = 0.6, end = 0.9),
    # Missing read counts
    date.seqres.data %>%
    mutate(`Has read n` = !is.na(`Read count`)) %>%
    ggplot() +
        aes(x = 1, fill = `Has read n`) +
        geom_bar(stat = "count") +
        geom_text(
            aes(label = after_stat(count)),
            stat = "count",
            position = position_stack(vjust = 0.5)
        ) +
        facet_grid("DateType") +
        theme(
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        scale_fill_viridis_d(begin = 0.6, end = 0.9),
    # Settings
    nrow = 1,
    align = "h",
    widths = c(8, 0.67, 0.67, 0.67),
    legend = "top"
)

# Building final plot
plots <- list(
    date.country.p,
    date.seqtech.p,
    date.seqres.p
)

ggarrange(plotlist = plots, ncol = 1)

message("Saving PDF report")
ggsave(
    snakemake@output$plot_pdf,
    width = snakemake@params$plot_width_in, height = 10 * length(plots),
    bg = "white"
)

message("Saving PNG report")
ggsave(
    snakemake@output$plot_png,
    width = snakemake@params$plot_width_in, height = 10 * length(plots),
    bg = "white"
)

message("Saving data tables")
write_csv(date.country.data, snakemake@output$country_timeline_table)
write_csv(date.seqtech.data, snakemake@output$tech_timeline_table)
write_csv(date.seqres.data, snakemake@output$seqres_timeline_table)
