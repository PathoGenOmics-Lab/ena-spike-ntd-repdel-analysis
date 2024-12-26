#!/usr/bin/env python

import csv
import logging
import sqlite3


if __name__ == "__main__":

    logging.basicConfig(
        level=logging.INFO,
        format=snakemake.config["PY_LOG_FMT"],
        filename=snakemake.log[0]
    )

    logging.info("Connecting to database")
    with open(snakemake.input.table) as f, sqlite3.connect(snakemake.output.database) as conn:
        reader = csv.DictReader(f, delimiter="\t")
        
        # Create table with text fields, using ROWID as PRIMARY KEY
        logging.info("Creating table ENARecords")
        columns = ",".join(f"{colname} TEXT" for colname in reader.fieldnames)
        conn.execute(f"CREATE TABLE ENARecords ({columns})")
        
        # Fill table
        logging.info("Filling table ENARecords")
        colnames = ",".join(reader.fieldnames)
        for record in reader:
            values = ",".join(f"'{record[colname]}'" for colname in reader.fieldnames)
            conn.execute(f"INSERT INTO ENARecords ({colnames}) VALUES ({values})")
