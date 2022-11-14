#!/usr/bin/env python3.9
import argparse
from pprint import pprint


def main():
    parser = argparse.ArgumentParser(
        prog='TreeChomp',
        description="Analyses your phylogenetic tree files to determine highest bootstrap\
                    for specified subtree criteria.")
    files_source = parser.add_mutually_exclusive_group()
    files_source.add_argument('-d', '--directory', default=".", dest="DIR",
                              help="Directory of .tre to analyze")
    files_source.add_argument('-l', '--list',
                              help="File containing list of files to analyze")
    files_source.add_argument('-f', '--file', nargs='+',
                              help=".tre files to analyze")

    continuation = parser.add_mutually_exclusive_group()
    continuation.add_argument('-r', action='store_true',
                              help="Replace previous output in the directory")
    continuation.add_argument('-c', action='store_true',
                              help="Attempt to continue previous output in the directory")

    parser.add_argument('-t', '--tolerance', default=0,
                        help="Absolute (1, 2) or relative (0.05, 0.1) tolerance of taxons outside of a set. Default is no tolerance.")

    parser.add_argument('criteria',
                        help="String list of criteria to apply on subtrees")
    parser.add_argument('sets',
                        help="String of definitions of taxons sets")

    parser.add_argument('-o', '--output',
                        help="Output CSV file")

    args = parser.parse_args()

    pprint(args)


if __name__ == "__main__":
    main()
