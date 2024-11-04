import sys
import logging
import time
import requests
from pathlib import Path

import pandas as pd


def download_file(url, path):
    r = requests.get(url)
    if r.status_code == 200:
        with open(path, "wb") as fw:
            fw.write(r.content)
        logging.debug(f"Downloaded {len(r.content)} bytes; sleeping {snakemake.params.sleep} s")
        time.sleep(snakemake.params.sleep)
    else:
        msg = f"Could not download {url} to {path} (status {r.status_code})"
        logging.error(msg)
        sys.exit(msg)


def format_url(text: str) -> str:
    return ("http://" + text) if not text.startswith("http://") else text


if __name__ == "__main__":

    logging.basicConfig(
        level=logging.INFO,
        format=snakemake.config["PY_LOG_FMT"],
        filename=snakemake.log[0]
    )
    
    # Read download index
    logging.info(f"Reading search records")
    run = pd.read_csv(snakemake.input.table)
    assert(len(run) == 1)
    
    # Download two FASTQ files
    urls = run.iloc[:, 0].fastq_ftp.split(";")
    assert(len(urls) == 1)
    url = format_url(urls[0])
    logging.info(f"Downloading {url} to {snakemake.output.fastq}")
    download_file(url, snakemake.output.fastq)
