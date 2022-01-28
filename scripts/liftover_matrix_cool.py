#!/usr/bin/env python3

import argparse

import bioframe as bf
import cooler
import h5py
import numpy as np
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser(description="Liftover-like script for .cool files.")
    cli.add_argument("--input-cool",
                     type=str,
                     required=True,
                     help="Path to input .cool (mcool URI syntax supported).")
    cli.add_argument("--output-cool",
                     type=str,
                     required=True,
                     help="Path to output .cool.")
    cli.add_argument("--chrom-sizes",
                     type=str,
                     required=True,
                     help="Path to .chrom.sizes file to use for producing the bin group for the output .cool file.")
    cli.add_argument("--chrom-ranges-bedpe",
                     type=str,
                     required=True,
                     help="Path to .bedpe with the genomic coordinate mappings to use during the liftover.")
    cli.add_argument("--assembly-name",
                     type=str,
                     required=True,
                     help="Assembly name used to generate .cool file(s)")
    return cli


def generate_bins_for_cooler(chroms, bin_size):
    bins = []
    for name, size in chroms.items():
        num_bins = (size + bin_size - 1) // bin_size
        bins.append(pd.DataFrame({"chrom": np.repeat(name, num_bins),
                                  "start": np.arange(0, size, bin_size, dtype=int),
                                  "end": np.clip(np.arange(0, size, bin_size, dtype=int) + bin_size, 0,
                                                 size)}))
    return pd.concat(bins, ignore_index=True)


def import_bedpe(path_to_bedpe):
    df = pd.read_table(path_to_bedpe,
                       names=["chrom1", "start1", "end1",
                              "chrom2", "start2", "end2"],
                       usecols=range(0, 6))
    return df.sort_values(by=["chrom1", "start1"])


def import_chrom_sizes(path_to_chrom_size):
    df = pd.read_table(path_to_chrom_size,
                       names=("chrom", "end"),
                       usecols=[0, 1])
    df["start"] = np.zeros_like(df["end"], dtype=int)
    return df[["chrom", "start", "end"]].reindex()


def import_chrom_sizes_from_cool(path_to_cool):
    chromsizes = cooler.Cooler(path_to_cool).chromsizes
    return pd.DataFrame({"chrom": chromsizes.index,
                         "start": np.zeros_like(len(chromsizes.index)),
                         "end": np.array(chromsizes)})


def read_pixels(path_to_cooler, chroms, starts, ends):
    """
    Return pixels for regions of interest only
    """
    assert len(chroms) == len(starts)
    assert len(starts) == len(ends)

    pixels = []
    c = cooler.Cooler(path_to_cooler)
    selector = c.matrix(balance=False, as_pixels=True)
    for chrom, start, end in zip(chroms, starts, ends):
        df = selector.fetch(f"{chrom}:{start}-{end}")
        df["chrom"] = chrom
        df["start"] = start
        pixels.append(df)

    return pd.concat(pixels, sort=True)[["chrom", "start", "bin1_id", "bin2_id", "count"]]


def fix_chrom_names(h5_handle):
    """
    Chromosome information under bins/ can be stored as enum/int.
    To simplify processing, replace chromosomes ids with their actual names
    """
    if isinstance(h5_handle["bins/chrom"].dtype, str):
        return h5_handle["bins/chrom"]

    names = {i: name for i, name in enumerate(h5_handle["chroms/name"])}
    return [names[i].decode() for i in h5_handle["bins/chrom"]]


def import_bins(path_to_cooler_file, bedpe_df, import_weights=True):
    """
    Here we are using h5py directly so that we can also import weights
    """
    toks = path_to_cooler_file.rsplit("::/resolutions/")
    root_group = "/"
    if len(toks) != 1:
        assert len(toks) == 2
        path_to_cooler_file = toks[0]
        root_group = "resolutions/" + toks[1].lstrip(":")

    with h5py.File(path_to_cooler_file, "r") as f:
        root = f[f"{root_group}/bins"]
        assert "chrom" in root.keys()
        assert "start" in root.keys()
        assert "end" in root.keys()

        if import_weights:
            # Import all datasets
            df1 = pd.DataFrame({k: root[k][:] for k in root.keys()})
        else:
            df1 = pd.DataFrame({k: root[k][:] for k in ["chrom", "start", "end"]})

        df1["chrom"] = fix_chrom_names(f[root_group])

    # Select bins overlapping regions of interest
    df2 = bf.overlap(df1,
                     bedpe_df,
                     return_index=True,
                     cols2=("chrom1", "start1", "end1"),
                     suffixes=("_1", "_2"))

    # NA entries in column index_2 indicate that no overlap was found for the corresponding row
    idx = df2.dropna(subset=["index_2"])["index_1"]

    return df1.loc[idx]


def compute_bin_offsets(bedpe_df, bin_size):
    offsets = (((bedpe_df["start2"] + bin_size - 1) // bin_size) -
               ((bedpe_df["start1"] + bin_size - 1) // bin_size)).astype(int)
    return pd.DataFrame({"chrom": bedpe_df["chrom1"],
                         "start": bedpe_df["start1"],
                         "offset": offsets})


def compute_offsets_in_bp(bin_offsets, bin_size):
    offsets = bin_offsets.copy()
    offsets["offset"] *= bin_size
    return offsets


def apply_offsets(df, label, offset_df):
    for (_, (chrom, start, offset)) in offset_df.iterrows():
        df.loc[(df["chrom"] == chrom) &
               (df["start"] == start), label] += offset

    return df


def join_bins_on_pixels_by_id(bins, pixels):
    bins1 = bins.copy()
    bins2 = bins.copy()
    # This will ensure that all columns will have unique (and predictable) identifiers
    bins1.columns = [f"{name}1" for name in bins1.columns]
    bins2.columns = [f"{name}2" for name in bins2.columns]

    # Rename bin_idN to binN_id
    bins1.columns = [name if name != "bin_id1" else "bin1_id" for name in bins1.columns]
    bins2.columns = [name if name != "bin_id2" else "bin2_id" for name in bins2.columns]

    # Join the two sets of bins
    pixels1 = pixels.merge(bins1,
                           how="left",
                           on="bin1_id",
                           suffixes=[None, None])
    return pixels1.merge(bins2,
                         how="left",
                         on="bin2_id",
                         suffixes=[None, None])


def join_bins_on_pixels_by_coords(bins, pixels):
    bins1 = bins.copy()
    bins2 = bins.copy()
    # This will ensure that all columns will have unique (and predictable) identifiers
    bins1.columns = [f"{name}1" for name in bins1.columns]
    bins2.columns = [f"{name}2" for name in bins2.columns]

    # Rename bin_idN to binN_id
    bins1.columns = [name if name != "bin_id1" else "bin1_id" for name in bins1.columns]
    bins2.columns = [name if name != "bin_id2" else "bin2_id" for name in bins2.columns]

    # Join the two sets of bins
    pixels1 = pixels.merge(bins1,
                           how="left",
                           on=["chrom1", "start1", "end1"],
                           suffixes=[None, None])
    return pixels1.merge(bins2,
                         how="left",
                         on=["chrom2", "start2", "end2"],
                         suffixes=[None, None])


def extract_weights_from_pixels_df(pixels):
    if len(pixels.columns) == 3:  # i.e. pixels df has no weights
        return None

    # Drop columns related to pixel information
    weights = pixels.drop(["chrom", "start", "count",
                           "chrom1", "start1", "end1",
                           "chrom2", "start2", "end2"],
                          axis="columns")

    # Split weights DF into weights1 and weights2. weights1 will contains weights referring to bin1_id while
    # weights2 will contain the same but for bin2_id
    weights1 = weights[[col for col in weights.columns if not col.endswith("2")]].drop("bin2_id", axis="columns")
    weights2 = weights[[col for col in weights.columns if not col.endswith("1")]].drop("bin1_id", axis="columns")
    weights1.columns = [name.removesuffix("1") for name in weights1.columns]
    weights2.columns = [name.removesuffix("2") for name in weights2.columns]

    # Remove duplicates and ensure that bin_ids can be used to uniquely identify rows
    weights1.drop_duplicates(inplace=True)
    weights2.drop_duplicates(inplace=True)
    assert len(weights1["bin1_id"].unique()) == len(weights1["bin1_id"])
    assert len(weights2["bin2_id"].unique()) == len(weights2["bin2_id"])

    # Use bin_ids as index
    weights1.index = weights1["bin1_id"]
    weights2.index = weights2["bin2_id"]

    weights1.drop("bin1_id", axis="columns", inplace=True)
    weights2.drop("bin2_id", axis="columns", inplace=True)

    # Concatenate weights and remove duplicate rows
    weights = pd.concat([weights1, weights2], axis="rows", sort=True)
    return weights[~weights.index.duplicated(keep="first")]


def clamp_pos(df, labels, chrom_sizes):
    if isinstance(labels, str):
        labels = [labels]
    for name, size in chrom_sizes.items():
        idx = df[df["chrom"] == name].index
        for label in labels:
            df.loc[idx, label] = np.clip(df.loc[idx, label], 0, size)
    return df


def copy_weights_attributes(input_cool_path, output_cool_path):
    toks = input_cool_path.rsplit("::/resolutions/")
    root_group1 = "/"
    if len(toks) != 1:
        assert len(toks) == 2
        input_cool_path = toks[0]
        root_group1 = "resolutions/" + toks[1].lstrip(":")

    toks = output_cool_path.rsplit("::/resolutions/")
    root_group2 = "/"
    if len(toks) != 1:
        assert len(toks) == 2
        output_cool_path = toks[0]
        root_group2 = "resolutions/" + toks[1].lstrip(":")

    with h5py.File(input_cool_path, "r") as f1, h5py.File(output_cool_path, "r+") as f2:
        if "weight" not in f1[f"{root_group1}/bins"]:
            return
        dset1 = f1[f"{root_group1}/bins/weight"]
        assert "weight" in f2[f"{root_group2}/bins"]
        dset2 = f2[f"{root_group2}/bins/weight"]

        for name, val in dset1.attrs.items():
            if name not in ["mad_max", "min_count", "min_nnz"]:
                dset2.attrs.create(name, val)


if __name__ == "__main__":
    args = make_cli().parse_args()

    # Parse/import inputs
    assembly_name = args.assembly_name
    bin_size = cooler.Cooler(args.input_cool).binsize

    chrom_sizes2 = {row["chrom"]: row["end"] for _, row in import_chrom_sizes(args.chrom_sizes).iterrows()}

    chrom_ranges = import_bedpe(args.chrom_ranges_bedpe)
    bin_offsets = compute_bin_offsets(chrom_ranges, bin_size)
    offsets = compute_offsets_in_bp(bin_offsets, bin_size)

    old_bins = import_bins(args.input_cool, chrom_ranges)
    old_bins["bin_id"] = old_bins.index

    old_pixels = read_pixels(args.input_cool,
                             chrom_ranges["chrom1"],
                             chrom_ranges["start1"],
                             chrom_ranges["end1"])

    # Join bin information to pixel DF (bin and pixel DF are both using the old gemonic coords)
    old_pixels = join_bins_on_pixels_by_id(old_bins, old_pixels)

    # Generate bins DF using the new coords
    new_bins = generate_bins_for_cooler(chrom_sizes2, bin_size)
    new_bins["bin_id"] = new_bins.index  # Storing a copy of the index simplifies downstream processing

    # Generate new pixels from old pixels
    new_pixels = old_pixels.copy()
    new_pixels = apply_offsets(new_pixels, "start1", offsets)
    new_pixels = apply_offsets(new_pixels, "start2", offsets)
    new_pixels["end1"] = new_pixels["start1"] + bin_size
    new_pixels["end2"] = new_pixels["start2"] + bin_size

    new_pixels = clamp_pos(new_pixels, ["end1", "end2"], chrom_sizes2)

    # Drop old bin_ids
    new_pixels.drop(["bin1_id", "bin2_id"], axis="columns", inplace=True)

    # Join new bin_ids on new pixels DF
    new_pixels = join_bins_on_pixels_by_coords(new_bins, new_pixels)

    # Extract new_weights from pixels DF.
    # Returns None if pixels DF does not have weights
    new_weights = extract_weights_from_pixels_df(new_pixels)

    # If present, drop weights columns from pixels DF
    new_pixels = new_pixels[["bin1_id", "bin2_id", "count"]]

    # Add weights to new_bins DF
    if new_weights is not None:
        for col in new_weights.columns:
            assert col not in new_bins
            new_bins[col] = np.nan

        new_bins.loc[new_weights.index, new_weights.columns] = new_weights

    # Drop unused column: bin_id is now new_bins' index
    new_bins.drop("bin_id", axis="columns", inplace=True)

    # Create new cooler file
    cooler.create_cooler(args.output_cool,
                         bins=new_bins,
                         pixels=new_pixels,
                         assembly=assembly_name)
    # ensure_sorted=True,
    # boundscheck=True,
    # dupcheck=True,
    # triucheck=True)

    copy_weights_attributes(args.input_cool, args.output_cool)
