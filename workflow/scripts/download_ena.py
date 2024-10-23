import sys
import time
import requests
from pathlib import Path

import pandas as pd


def download_file(url, path):
    r = requests.get(url)
    if r.status_code == 200:
        with open(path, "wb") as fw:
            fw.write(r.content)
        time.sleep(snakemake.params.sleep)
    else:
        sys.exit(f"Could not download {url} to {path} (status {r.status_code})")


def format_url(text: str) -> str:
    return ("http://" + text) if not text.startswith("http://") else text


if __name__ == "__main__":
    
    # Read download index
    run = pd.read_csv(snakemake.input.table)
    print(f"Read {len(run)} search records", flush=True)

    # Prepare file download
    print("Downloading data", flush=True)
    folder = Path(snakemake.output.folder)
    folder.mkdir(parents=True)
    
    # Perform file download
    for i, row in run.iterrows():
        for url in row.fastq_ftp.split(";"):
            url = format_url(url)
            file_name = Path(url).name
            print(f"Downloading {url} to {folder / file_name}", flush=True)
            download_file(url, folder / file_name)
