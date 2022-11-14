#!/usr/bin/env python3.9
import argparse
from pprint import pprint
from os import listdir
from os.path import isfile, join


def main():
    args = parse_args()
    pprint(args)

    io = IOHandler(directory=args.directory, files=args.files, flist=args.list, output=args.output)

    # while f := io.read_file():
    #     print(f)


class IOHandler:
    default_csv = "tc-output.csv"

    @staticmethod
    def extension_is_csv(filename: str):
        return len(filename) > 4 and filename[-4:].lower() == '.csv'

    def __init__(self, directory=None, files=None, flist=None, output=None):
        # get a list of valide files to analyze
        self.tree_files = []
        if directory:
            directory = directory[0]
            # use all files in the set directory
            self.directory = directory
            self.tree_files = [f for f in listdir(directory) if isfile(join(directory, f))]
        elif files:
            # use files listed as arguments
            self.tree_files = [f for f in files if isfile(f)]
        elif flist:
            flist = flist[0]
            # use list-file to get files listed in it
            if isfile(flist):
                with open(flist) as file:
                    potential_files = [line.rstrip() for line in file]
                    self.tree_files = [f for f in potential_files if isfile(f)]

        self.tree_files = list(filter(lambda x: not self.extension_is_csv(x), self.tree_files))
        self.tree_files.sort()

        self.current_file_index = 0
        self.file_count = len(self.tree_files)

        if output:
            output = output[0]

            if self.extension_is_csv(output):
                self.output = output
            else:
                self.output = output + ".csv"
        else:
            self.output = self.default_csv

    def read_file(self, n=None):
        if n and n < self.file_count:
            with open(self.tree_files[n]) as file:
                return file.readlines()[0]
        elif self.current_file_index < self.file_count:
            # just read next
            with open(self.tree_files[self.current_file_index]) as file:
                self.current_file_index += 1
                return file.readlines()[0]
        else:
            return None

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
