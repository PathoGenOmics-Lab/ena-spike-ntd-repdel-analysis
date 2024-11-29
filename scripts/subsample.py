import sys
import random
import argparse
import logging


def count_lines(path: str) -> int:
    with open(path) as f:
        return sum(1 for _ in f)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("--subsample", type=int, default=0)
    parser.add_argument("--subsample-seed", type=int, default=7291)
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    logging.info("Counting input records (excluding header)")
    n_records = count_lines(args.input) - 1

    if n_records < args.subsample:
        logging.error(f"The number of records (n={n_records}) is lower than the requested subsample size (n={args.subsample})")
        sys.exit(1)

    logging.info("Calculating subsample")
    random.seed(args.subsample_seed)
    indices = list(range(n_records))
    random.shuffle(indices)
    indices = indices[:args.subsample]
    
    logging.info("Writing subsample")
    with open(args.input) as f, open(args.output, "w") as fw:
        # Write header
        fw.write(f.readline())
        # Write rows from subsample
        for i, line in enumerate(f):
            if i in indices:
                fw.write(line)
