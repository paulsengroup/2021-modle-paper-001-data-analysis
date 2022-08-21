#!/usr/bin/env python3

import argparse
import bioframe as bf
import sys


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument("bins")
    cli.add_argument("genome_assembly")

    return cli


if __name__ == "__main__":
    args = vars(make_cli().parse_args())

    df = bf.frac_gc(bf.read_table(args["bins"], schema="bed3"),
                    bf.load_fasta(args["genome_assembly"]))

    df.to_csv(sys.stdout, index=False, header=False, sep="\t")
