#!/usr/bin/env python

import csv
import logging
import argparse
import sqlite3


if __name__ == "__main__":

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input TSV")
    parser.add_argument("output", help="Output SQLite database")
    args = parser.parse_args()

    logging.info("Connecting to database")
    with open(args.input) as f, sqlite3.connect(args.output, isolation_level=None) as conn:
        cursor = conn.execute("PRAGMA journal_mode=WAL")
        if cursor:
            logging.info(f"Database set to journal mode = {'/'.join(cursor.fetchone())}")
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
