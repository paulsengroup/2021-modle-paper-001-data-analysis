#!/usr/bin/env python3

import os
import re
import shutil
import argparse
import cooler


def make_cli():
    cli = argparse.ArgumentParser()
    cli.add_argument("coolers",
                     nargs="+",
                     type=str,
                     help="One or more .(m)cool files whose chromosome names should be normalized.\n"
                          "Chromosomes are renamed such that:\n"
                          "- 1 -> chr1\n"
                          "- chrom1 -> chr1\n"
                          "- chr1 -> chr1")
    cli.add_argument("--inplace",
                     default=False,
                     action="store_true",
                     help="Rename chromosomes inplace.\n"
                          "By default a copy of the input file(s) is made and the new files are modified.\n"
                          "New files are stored in the same location as input files (and have suffix .new)")
    return cli


if __name__ == "__main__":
    args = vars(make_cli().parse_args())
    pattern = re.compile(r"^chrom|chr", re.IGNORECASE)

    for input_mcool in args["coolers"]:
        output_mcool = os.path.basename(input_mcool.strip()) + ".new"
        shutil.copyfile(input_mcool, output_mcool)

        for path in cooler.fileops.list_coolers(input_mcool):
            c = cooler.Cooler(f"{output_mcool}::{path}")
            mappings = {chrom: "chr" + pattern.sub("", chrom, 1) for chrom in c.chromnames}
            cooler.rename_chroms(c, mappings)
