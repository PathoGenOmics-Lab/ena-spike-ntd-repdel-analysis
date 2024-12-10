import argparse
import logging

import pandas as pd
from typing import Iterator, List


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("--seed", type=int, default=7291)
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    logging.info("Shuffling")
    pd.read_csv(args.input, sep="\t") \
        .sample(frac=1, random_state=args.seed) \
        .to_csv(args.output, index=False)
