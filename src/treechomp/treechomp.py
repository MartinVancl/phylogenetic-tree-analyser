#!/usr/bin/env python3.9
import argparse
from pprint import pprint
from os import listdir
from os.path import isfile, join


def main():
    args = parse_args()
    pprint(args)

    io = IOHandler(directory=args.directory, files=args.files, flist=args.list, output=args.output)
    print("IO TREE_FILES:")
    pprint(io.tree_files)
    pprint(io.output)
    print(f"{args.criteria=}")
    print(f"{args.sets=}")


class IOHandler:
    default_csv = "tc-output.csv"

    def __init__(self, directory=None, files=None, flist=None, output=None):
        #
        self.tree_files = []
        if directory:
            # use all files in the set directory
            self.directory = directory
            self.tree_files = [f for f in listdir(directory) if isfile(join(directory, f))]
        elif files:
            # use files listed as arguments
            self.tree_files = [f for f in files if isfile(f)]
        elif flist:
            # use list-file to get files listed in it
            if isfile(flist):
                with open(flist) as file:
                    potential_files = [line.rstrip() for line in file]
                    self.tree_files = [f for f in potential_files if isfile(f)]

        if output:
            output = output[0]

            if len(output) > 4 and output[-4, 0] == ".csv":
                self.output = output
            else:
                self.output = output + ".csv"
        else:
            self.output = self.default_csv

    class FileLoader:
        pass

    class CSVExport:
        pass

    pass


class PTree:
    class Edge:
        pass

    class Node:
        pass

    class Taxon(Node):
        pass


class TreeFile:
    def __init__(self, filename):
        # TODO takes FASTA format and builds a tree
        pass


def parse_args():
    parser = argparse.ArgumentParser(
        prog='TreeChomp',
        description="Analyses your phylogenetic tree files to determine highest bootstrap\
                        for specified subtree criteria.")
    files_source = parser.add_mutually_exclusive_group()
    files_source.add_argument('-d', '--directory', nargs=1,
                              help="Directory of tree files to analyze")
    files_source.add_argument('-l', '--list', nargs=1,
                              help="Files containing list of files to analyze")
    files_source.add_argument('-f', '--files', nargs='*',
                              help=".tre files to analyze")

    continuation = parser.add_mutually_exclusive_group()
    continuation.add_argument('-n', action='store_true', default=True,
                              help="Start new run and replace previous CSV output in the directory")
    continuation.add_argument('-a', action='store_true',
                              help="Attempt to append to previous CSV output in the directory")

    parser.add_argument('-t', '--tolerance', default=0,
                        help="Absolute (1, 2) or relative (0.05, 0.1) tolerance \
                            of taxons outside of a set. Default is no tolerance.")

    parser.add_argument('-o', '--output', nargs=1,
                        help="Output CSV file")
    parser.add_argument('-v', action='store_true',
                        help="Run verbose")

    parser.add_argument('criteria', nargs=1,
                        help="String list of criteria to apply on subtrees")
    parser.add_argument('sets', nargs=1,
                        help="String of definitions of taxons sets")

    return parser.parse_args()


if __name__ == "__main__":
    main()
