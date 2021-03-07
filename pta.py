#!/usr/bin/env python3

from sys import argv
from os import path, listdir
from os.path import isfile, join

class argParser:
    def __init__(self, verbose=False):
        self.fullArgsString = argv[0]
        self.args = argv[1:]
        self.argsLength = len(self.args)

        if verbose:
            print("args=", end='')
            print(self.args)

    def isFlagSet(self, flag):
        return f"-{flag}" in self.args

    def getFlagPosition(self, flag):
        if self.isFlagSet(flag):
            return self.args.index(f"-{flag}")
        else:
            return None

    def getFlagValue(self, flag):
        fv = self.isFlagValued(flag, True)
        if fv[0]:
            return self.args[fv[1]]
        else:
            return None
            
    def isFlagValued(self, flag, giveIndex=False):
        i = self.getFlagPosition(flag)
        if type(i) == int:
            if i + 1 < self.argsLength and self.args[i + 1][0] != '-':
                # there is an argument after the flag and it is not another flag
                return (True, i + 1) if giveIndex else True
            else:
                # flag is not valued with argument
                return (False, None) if giveIndex else False
        else:
            #flag is not set
            return (False, None) if giveIndex else False

class treeFilesManager:
    def __init__(self, args):
        self.args = args
        if self.args.isFlagValued('d'):
            # working directory is set by -d ...
            wd = self.args.getFlagValue('d')
            if path.exists(wd):
                self.wd = wd
            else:
                exit(f"Working directory '{wd}' set by -d does not exist.")
        else:
            self.wd = '.'

    def fileExists(self, filename):
        return isfile(join(self.wd, filename))

    def loadFromFile(self, filename):
        pass
    def loadFromDirectory(self, directory):
        pass
    def filterExisting(self):
        pass

def main():

    args = argParser()

    treeFiles = treeFilesManager(args)

    print("exists", treeFiles.fileExists("Takayama-helix-ANV-100120_homologs.fasta.tre"))



def noexec():
    # tis an o'code that i dissect and transform into something more OOP

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