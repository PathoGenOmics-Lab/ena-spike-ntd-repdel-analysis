import logging
from pathlib import Path
from typing import List, Mapping

import pandas as pd


# input example:
# CHROM	        REF	        POS	    ALT	DP	    ALT_DP	ALT_RV	ALT_FREQ	ALT_QUAL	GENE	                                        HGVS_P
# NC_045512.2	ATACATG	    21764	A	1160	512	    0	    0.441379	20	        S,ORF3a,E,M,ORF1ab,ORF1ab,ORF1ab,ORF1ab,ORF1ab	p.His69_Val70del,,,,,,,,
# NC_045512.2	GGTGTTTATT	21986	G	5142	2843	0	    0.552898	20	        S,ORF3a,E,M,ORF1ab,ORF1ab,ORF1ab,ORF1ab,ORF1ab	p.Gly142_Tyr145delinsAsp,,,,,,,,


def pass_exclude(df: pd.DataFrame, markers: List[Mapping[str, str]], freq_threshold: float) -> bool:
    checks = []
    for marker in markers:
        gene, expression = marker["gene"], marker["expression"]
        selection = df[(df.GENE == gene) & (df.HGVS_P == expression) & (df.ALT_FREQ <= freq_threshold)]
        logging.info(f"Found {len(selection)} variants for {gene}:{expression} with frequency <= {freq_threshold}")
        checks.append(len(selection) == 0)
    pass_checks = all(checks)
    logging.info(f"Pass exclude: {pass_checks}")
    return pass_checks


def pass_include(df: pd.DataFrame, markers: List[Mapping[str, str]], freq_threshold: float) -> bool:
    checks = []
    for marker in markers:
        gene, expression = marker["gene"], marker["expression"]
        selection = df[(df.GENE == gene) & (df.HGVS_P == expression) & (df.ALT_FREQ >= freq_threshold)]
        logging.info(f"Found {len(selection)} variants for {gene}:{expression} with frequency >= {freq_threshold}")
        checks.append(len(selection) >= 1)
    pass_checks = all(checks)
    logging.info(f"Pass include: {pass_checks}")
    return pass_checks


def is_haplotype(df: pd.DataFrame, marker_bunch: Mapping[str, List[Mapping[str, str]]], min_include_freq: float, max_exclude_freq: float) -> bool:
    return pass_exclude(df, marker_bunch["exclude_hgvs_p"], max_exclude_freq) and pass_include(df, marker_bunch["include_hgvs_p"], min_include_freq)


def read_variants(path: str, columns: List[str]) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df.columns = columns
    annotation = df.apply(lambda row: list(zip(row["GENE"].split(","), row["HGVS_P"].split(","))), axis=1) \
        .explode() \
        .apply(lambda x: pd.Series(x, index=["GENE", "HGVS_P"])) \
        .groupby(level=0).ffill()
    return  df.drop(["GENE", "HGVS_P"], axis=1).join(annotation)


def extract_sample_fields(path: str) -> Mapping[str, str]:
    path = Path(path)
    parts = path.parts
    print(parts)
    fields = {
        "study": parts[3],
        "sample": parts[4],
        "platform": parts[5],
        "run": parts[6],
        "haplotype": path.stem
    }
    fields["layout"], fields["nfastq"], fields["strategy"] = parts[7].split("_")
    return fields


if __name__ == "__main__":

    logging.basicConfig(
        level=logging.INFO,
        format=snakemake.config["PY_LOG_FMT"],
        filename=snakemake.log[0]
    )

    # Read input params
    unique_markers = set(f"{marker['gene']}:{marker['expression']}" for key in snakemake.params.markers.keys() for marker in snakemake.params.markers[key])
    logging.info(f"Unique markers: {unique_markers}")
    min_include_freq = float(snakemake.wildcards.inclpct) / 100
    max_exclude_freq = float(snakemake.wildcards.exclpct) / 100
    logging.info(f"Frequency thresholds: min include = {min_include_freq}; max exclude = {max_exclude_freq}")

    i = 0
    results = pd.DataFrame(columns=("study", "sample", "platform", "run", "layout", "nfastq", "strategy", "haplotype"))
    for variants_path in snakemake.input:
        logging.info(f"Processing {variants_path}")
        variants = read_variants(variants_path, snakemake.params.columns)
        if is_haplotype(variants, snakemake.params.markers, min_include_freq, max_exclude_freq):
            logging.debug(f"Filling positive results for {variants_path}")
            # Add sample identifiers
            results.loc[i] = extract_sample_fields(variants_path)
            # Add marker data from variant calling
            for _, row in variants.iterrows():
                formatted_marker = row["GENE"] + ":" + row["HGVS_P"]
                logging.debug(f"Adding result column for marker {formatted_marker}")
                if formatted_marker in unique_markers:
                    results.loc[i, formatted_marker] = f"{row['POS']}{row['REF']}>{row['ALT']}," + \
                        ",".join(str(it) for it in (row["DP"], row["ALT_DP"], row["ALT_RV"], row["ALT_FREQ"], row["ALT_QUAL"]))
            i += 1
    
    logging.info(f"Writing {i} records")
    results.to_csv(snakemake.output.table, index=False)
