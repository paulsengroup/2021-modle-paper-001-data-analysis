import argparse

from PolymerCpp.helpers import getCppSAWLC


def check_uint(value):
    x = int(value)
    if x <= 0:
        raise argparse.ArgumentTypeError(f"{x} is not a valid positive integer value")
    return x


def check_bounds_fp(value, lb=0.0, ub=1.0):
    x = float(value)
    if x < lb or x > ub:
        raise argparse.ArgumentTypeError(f"{x} is not between {lb}, {ub}")
    return x


def make_cli():
    cli = argparse.ArgumentParser()
    cli.add_argument("--number-of-monomers",
                     required=True,
                     type=check_uint)
    cli.add_argument("--persistence-length",
                     default=5,
                     type=check_uint)
    cli.add_argument("--link-diameter",
                     default=1.0,
                     type=check_bounds_fp)
    cli.add_argument("--seed",
                     default=883107317964742334,
                     type=check_uint)
    return cli


if __name__ == "__main__":
    args = make_cli().parse_args()

    assert args.number_of_monomers > 1
    coords = getCppSAWLC(args.number_of_monomers - 1,
                         args.persistence_length,
                         args.link_diameter,
                         args.seed)

    print(f"{args.number_of_monomers}\n")
    for x, y, z in coords:
        print(f"A\t{x:.18g}\t{y:.18g}\t{z:.18g}")
