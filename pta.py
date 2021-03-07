#!/usr/bin/env python3

from sys import argv
from os import path, listdir
from os.path import isfile, join

def main():
    args = argv[1:]
    print("args=", end='')
    print(args)
    argsLen = len(args)


    # check input file argument -f
    treesFile = False
    if '-f' in args:
        i = args.index('-f')
        if i + 1 < argsLen:
            fpath = args[i+1]
            if path.exists(fpath):
                treesFile = fpath
                # print(f"File exists: '{treesFile}'")
            else:
                print(f"Input file '{fpath}' does not exist.")
                exit()

    # check output file argument -o
    if '-o' in args:
        i = args.index('-o')
        if i + 1 < argsLen:
            fpath = args[i+1]
            if path.exists(fpath):
                if '-a' in args and '-r' not in args:
                    print("File exists and will be appended")
                elif '-r' in args and '-a' not in args:
                    print("File exists and will be replaced")
                else:
                    print(f"File {fpath} already exists. \nUse either -a flag to append results or -r to replace results (WARNING! this deletes content of '{fpath}').")


    folderpath = '.'
    if not treesFile:
        treesFile = [f for f in listdir(folderpath) if isfile(join(folderpath, f)) and f[-9:] == 'fasta.tre']
    else:
        with open(fpath, 'r') as file:
            treesFile = []
            for line in file.readlines():
                
                # if line[-9:] == 'fasta.tre' or line[-10:-1] == 'fasta.tre':
                treesFile += [line.rstrip()]
            count = len(treesFile)
            treesFile = [f for f in treesFile if isfile(f) and f[-9:] == 'fasta.tre']
            deleted = count - len(treesFile)
            if deleted > 0:
                print(f"{deleted} file(s) in '{fpath}' don't exist and are ignored")

    print(treesFile)








if __name__ == "__main__":
    main()