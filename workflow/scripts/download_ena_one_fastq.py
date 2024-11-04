import sys
import logging
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

import pandas as pd


def download_file(session, url, path):
    r = session.get(url)
    if r.status_code == 200:
        with open(path, "wb") as fw:
            fw.write(r.content)
        logging.debug(f"Downloaded {len(r.content)} bytes; sleeping {snakemake.params.sleep} s")
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

    session = requests.Session()
    retry = Retry(
        connect=snakemake.params.retries,
        backoff_factor=snakemake.params.sleep
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("http://", adapter)
    
    # Read download index
    logging.info("Reading search records")
    run = pd.read_csv(snakemake.input.table)
    assert(len(run) == 1)
    
    # Download two FASTQ files
    logging.info("Formatting URL")
    urls = run.iloc[0].fastq_ftp.split(";")
    assert(len(urls) == 1)
    url = format_url(urls[0])
    logging.info(f"Downloading {url} to {snakemake.output.fastq}")
    download_file(session, url, snakemake.output.fastq)
