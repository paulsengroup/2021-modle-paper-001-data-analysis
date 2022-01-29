#!/usr/bin/env python3

import argparse
import itertools
from pathlib import Path

import cooler
import numpy as np
import pandas as pd
from PIL import Image, ImageOps

def make_cli():
    cli = argparse.ArgumentParser()
    cli.add_argument("--input-pictures",
                     nargs="+",
                     type=str,
                     required=True,
                     help="Path to one or more pictures to be processed")
    cli.add_argument("--chrom-ranges-bed",
                     type=str,
                     required=True,
                     help="Path to a BED file with the genomic coordinates of the regions shown in the input picture(s).")
    cli.add_argument("--chrom-sizes",
                     type=str,
                     required=True,
                     help="Path to .chrom.size file to use when generating .cool files.")
    cli.add_argument("--bin-size",
                     type=int,
                     required=True,
                     help="Bin size used to generate .cool file(s).")
    cli.add_argument("--assembly-name",
                     type=str,
                     required=True,
                     help="Assembly name used to generate .cool file(s)")
    cli.add_argument("--output-dir",
                     required=True,
                     type=str)
    return cli


def import_image(path):
    img = ImageOps.invert(Image.open(path).convert('L'))
    return np.array(img, dtype=int)


def resize_image(pixels, new_size):
    assert pixels.shape[0] == pixels.shape[1]
    return np.array(Image.fromarray(pixels).resize([new_size, new_size],
                                                   resample=Image.LANCZOS,
                                                   reducing_gap=3), dtype=int)


def generate_bins_for_cooler(chroms, bin_size):
    bins = []
    for name, size in chroms.items():
        num_bins = (size + bin_size - 1) // bin_size
        bins.append(pd.DataFrame({"chrom": np.repeat(name, num_bins),
                                  "start": np.arange(0, size, bin_size, dtype=int),
                                  "end": np.clip(np.arange(0, size, bin_size, dtype=int) + bin_size, 0,
                                                 size)}))
    return pd.concat(bins, ignore_index=True)


def get_bin_ids(bins, chrom_name, start, end):
    return bins[(bins["chrom"] == chrom_name) &
                (bins["start"] >= start) &
                (bins["start"] < end)].index.to_numpy()


def generate_cooler_pixels_from_img(img, chrom_name, start, end, bins):
    region_size = end - start
    assert region_size > 0

    bin_size = int(bins.head(1)["end"] - bins.head(1)["start"])
    num_bins_region = int(round(region_size / bin_size))

    bin1_ids = np.array(list(itertools.product(get_bin_ids(bins, chrom_name, start, end), repeat=2)))
    bin2_ids = bin1_ids[:, 1]
    bin1_ids = bin1_ids[:, 0]
    pixels = resize_image(np.array(img, dtype=np.uint8), num_bins_region).flatten()
    pixels -= np.min(pixels)

    print("bin1", bin1_ids.shape)
    print("bin2", bin2_ids.shape)
    print("pixels", pixels.shape)
    print("expected_bins", num_bins_region)
    assert bin1_ids.shape == pixels.shape
    assert bin2_ids.shape == pixels.shape

    pixels = pd.DataFrame({"bin1_id": bin1_ids,
                           "bin2_id": bin2_ids,
                           "count": pixels})
    pixels = pixels[pixels["count"] != 0]
    return pixels[pixels["bin1_id"] <= pixels["bin2_id"]]


def make_cooler(path_to_chrom_sizes, path_to_chrom_ranges, path_to_images, path_to_output_cool, assembly_name,
                bin_size):
    chrom_sizes = pd.read_table(path_to_chrom_sizes,
                                header=None,
                                names=("chrom", "size"))
    chrom_ranges = pd.read_table(path_to_chrom_ranges,
                                 names=("chrom", "start", "end"),
                                 usecols=range(0, 3))

    assert len(chrom_ranges) == len(path_to_images)

    print("Generating bins...")
    bins = generate_bins_for_cooler({name: size for (_, (name, size)) in chrom_sizes.iterrows()}, bin_size)
    print("Importing images...")
    imgs = [import_image(path) for path in path_to_images]
    print("Processing images...")
    pixels = [generate_cooler_pixels_from_img(img, name, start, end, bins) for img, name, start, end in
              zip(imgs, chrom_ranges["chrom"], chrom_ranges["start"], chrom_ranges["end"])]

    print("Aggregating contacts...")
    pixels = pd.concat(pixels, sort=True)

    print("Generating cooler...")
    cooler.create_cooler(path_to_output_cool, bins, pixels, assembly=assembly_name)


if __name__ == "__main__":
    args = make_cli().parse_args()
    out_dir = args.output_dir
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    make_cooler(args.chrom_sizes,
                args.chrom_ranges_bed,
                args.input_pictures,
                f"{out_dir}/{args.assembly_name}_{args.bin_size}.cool",
                args.assembly_name,
                args.bin_size)
