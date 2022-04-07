#!/usr/bin/env python3

import argparse
from sys import stdout

import pandas as pd


def make_cli():
    class SplitCSVArgs(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            setattr(namespace, self.dest, values.split(","))

    cli = argparse.ArgumentParser()
    cli.add_argument("benchmark_report_tsv",
                     type=str,
                     help="Path to benchmark report in TSV format.")
    cli.add_argument("--groupby-keys",
                     required=True,
                     type=str,
                     action=SplitCSVArgs,
                     help="Comma separated list of column name(s) to use for aggregation.")
    cli.add_argument("--skip-columns",
                     type=str,
                     action=SplitCSVArgs,
                     default=["id", "date", "hostname"],
                     help="Comma separated list of columns to exclude from the summary.")
    cli.add_argument("--method",
                     type=str,
                     default="median",
                     choices={"mean", "median"},
                     help="Aggregation method.")
    return cli


if __name__ == "__main__":
    args = make_cli().parse_args()

    aggr_method = args.method
    report = pd.read_table(args.benchmark_report_tsv)

    if args.skip_columns and isinstance(args.skip_columns, list) and len(args.skip_columns) != 0:
        report.drop(args.skip_columns, axis="columns", errors="ignore", inplace=True)

    if "max_resident_mem_kb_gpu" in report:
        report["max_resident_mem_kb"] += report["max_resident_mem_kb_gpu"]
        report.drop("max_resident_mem_kb_gpu", axis="columns")


    def aggregate(data, keys, method):
        if method == "median":
            return data.groupby(by=keys).median()
        if method == "std":
            return data.groupby(by=keys).std()

        assert method == "mean"
        return data.groupby(by=keys).mean()


    summary = pd.merge(aggregate(report, args.groupby_keys, aggr_method),
                       aggregate(report, args.groupby_keys, "std"),
                       left_index=True,
                       right_index=True,
                       suffixes=(f"_{aggr_method}", "_std"),
                       sort=True).reset_index()
    cols_to_drop = {k.removesuffix("_std") for k in summary.columns if
                    k.endswith("_std") and len(summary[k].unique()) == 1}

    summary.drop([f"{k}_std" for k in cols_to_drop], axis="columns", inplace=True)

    for col in cols_to_drop:
        summary.rename({f"{col}_{aggr_method}": col}, inplace=True)

    for col in summary.select_dtypes(include="float64"):
        if summary[col].apply(lambda n: n.is_integer()).all():
            summary[col] = summary[col].astype(int)

    summary.to_csv(stdout, sep="\t", index=False, float_format="%.2f")
