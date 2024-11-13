SEARCH_DF_COLS = (
    "study_accession", "sample_accession", "instrument_platform",
    "run_accession", "library_layout", "fastq_ftp", "library_strategy"
)


def count_fastq(row: dict) -> int:
    return row["fastq_ftp"].count(";") + 1


def build_groups(wildcards, table, columns) -> list:
    delimiter = "\t" if table.endswith(".tsv") else ","
    with open(table) as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        groups = {
            tuple(row[col] if col != "fastq_ftp" else count_fastq(row) for col in columns) \
            for row in reader
        }
    return sorted(groups)


def build_groups_filtering(wildcards, table, columns, **kwargs) -> list:
    groups = set()
    for group in build_groups(wildcards, table, columns):
        filters = []
        for column, values in kwargs.items():
            if type(values) is not list or type(value) is not tuple:
                values = [values]
            index = columns.index(column)
            filters.append(group[index] in values)
        if all(filters):
            groups.add(group)
    return groups


def build_search_targets(wildcards, template: str, columns=SEARCH_DF_COLS) -> list:
    return sorted(
        template.format(*groups) \
        for groups in build_groups(
            wildcards,
            checkpoints.filter_search_ena.get(**wildcards).output.table,
            columns
        )
    )


def build_search_targets_filtering(wildcards, template: str, columns=SEARCH_DF_COLS, **kwargs) -> list:
    columns = tuple(list(columns) + [column for column in kwargs.keys()])
    return sorted(
        template.format(*groups) \
        for groups in build_groups_filtering(
            wildcards,
            checkpoints.filter_search_ena.get(**wildcards).output.table,
            columns,
            **kwargs
        )
    )


def build_pangolin_targets(wildcards, template: str, columns=SEARCH_DF_COLS) -> list:
    return sorted(
        template.format(*groups) \
        for groups in build_groups(
            wildcards,
            checkpoints.filter_search_ena_with_pangolin.get(**wildcards).output.search_table,
            columns
        )
    )


def build_pangolin_targets_filtering(wildcards, template: str, columns=SEARCH_DF_COLS, **kwargs) -> list:
    columns = tuple(list(columns) + [column for column in kwargs.keys()])
    return sorted(
        template.format(*groups) \
        for groups in build_groups_filtering(
            wildcards,
            checkpoints.filter_search_ena_with_pangolin.get(**wildcards).output.search_table,
            columns,
            **kwargs
        )
    )


def build_snpsift_hgvs_p_filter(wildcards):
    items = set()
    for marker_class in ("include_hgvs_p", "exclude_hgvs_p"):
        for marker in config["HAPLOTYPES"][w.haplotype].get(marker_class, []):
            gene, expression = marker["gene"], marker["expression"]
            items.add(f"(ANN[*].GENE = '{gene}' & ANN[*].HGVS_P = '{expression}')")
    if len(items) > 0:
        return f" & ({' | '.join(items)})"
    else:
        return ""
