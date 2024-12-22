#!/usr/bin/env python


if __name__ == "__main__":
    with open(snakemake.output, "w") as fw:
        for path in snakemake.input:
            with open(path) as f:
                for line in f:
                    fw.write(line)
