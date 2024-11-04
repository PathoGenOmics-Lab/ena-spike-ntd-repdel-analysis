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
        format=snakemake.config.PY_LOG_FMT,
        filename=snakemake.log
    )
    
    # Read download index
    run = pd.read_csv(snakemake.input.table)
    logging.info(f"Read {len(run)} search records")

    # Prepare file download
    logging.info("Downloading data")
    folder = Path(snakemake.output.folder)
    folder.mkdir(parents=True)
    
    # Perform file download
    for i, row in run.iterrows():
        for url in row.fastq_ftp.split(";"):
            url = format_url(url)
            file_name = Path(url).name
            logging.info(f"Downloading {url} to {folder / file_name}")
            download_file(url, folder / file_name)
