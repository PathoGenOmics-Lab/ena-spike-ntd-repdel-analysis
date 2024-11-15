import logging
import re
import pandas as pd


COLUMNS = [
        "study_accession", "sample_accession", "instrument_platform",
        "run_accession", "library_layout", "library_strategy"
    ]

CONSENSUS_PATTERN = re.compile(snakemake.params.consensus_pattern)
COVERAGE_PATTERN = re.compile(snakemake.params.coverage_sample_pattern)


if __name__ == "__main__":

    logging.basicConfig(
        level=logging.INFO,
        format=snakemake.config["PY_LOG_FMT"],
        filename=snakemake.log[0]
    )

    # Read pangolin assignment
    logging.info("Reading pangolin assignment")
    pangolin = pd.read_csv(snakemake.input.pangolin_table)
    logging.info(f"Read {len(pangolin)} records")

    # Extract sample metadata from record ID
    logging.info("Extracting metadata from taxons")
    if len(pangolin) != 0:
        pangolin[COLUMNS] = pangolin["taxon"].apply(lambda record_id: CONSENSUS_PATTERN.match(record_id).groups()).tolist()
        # Keep only sample metadata
        pangolin = pangolin.drop([col for col in COLUMNS if col not in COLUMNS], axis="columns")
    
    # Read coverage
    logging.info("Reading coverage data")
    coverage = pd.read_csv(snakemake.input.coverage_table)
    logging.info(f"Read {len(coverage)} records")

    # Extract sample metadata from sample field
    logging.info("Extracting metadata from sample field")
    if len(coverage) != 0:
        coverage[COLUMNS] = coverage["sample"].apply(lambda sample: COVERAGE_PATTERN.match(sample).groups()).tolist()
        # Keep only sample metadata
        coverage = coverage.drop([col for col in COLUMNS if col not in COLUMNS], axis="columns")

    # Read filtered search results
    logging.info("Reading filtered search results")
    search = pd.read_csv(snakemake.input.search_table, sep="\t")
    
    # Keep only search results in the pangolin table and write
    logging.info("Filtering and writing")
    search \
        .merge(pangolin, how="right", on=COLUMNS) \
        .merge(coverage, how="right", on=COLUMNS) \
        .to_csv(snakemake.output.search_table, sep="\t", index=False)
