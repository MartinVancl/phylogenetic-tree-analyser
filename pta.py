#!/usr/bin/env python3

DEFAULT_OUTPUT_FILE = 'pta_output.csv' # for unspecified output file
TREE_FILE_EXTENSION = '.fasta.tre' # used when a whole directory of trees is set
RUN_QUIETLY = False # no prints will be executed

if RUN_QUIETLY:
    # this overrides the built-in print function to silence the whole program
    def overriden_print(*args):
        pass
    print = overriden_print

from sys import argv
from os import path, listdir
from os.path import isfile, join
from re import split, finditer, sub

class argParser:
    def __init__(self, verbose=False):
        # load arguments from system interface
        # set basic variables within the instance
        self.fullArgsString = argv[0]
        self.args = argv[1:]
        self.argsLength = len(self.args)

        if verbose:
            print("args=", end='')
            print(self.args)

    def isFlagSet(self, flag):
        # boolean representation of whether a flag is present in arguments
        return f"-{flag}" in self.args

    def getFlagPosition(self, flag):
        # returns position of a flag as index, nonexistent returns None
        if self.isFlagSet(flag):
            # flag exists, find it's index
            return self.args.index(f"-{flag}")
        else:
            return None

    def getFlagValue(self, flag):
        # for a flag followed by a value, gives the value or None for non-existent flag or value
        fv = self.isFlagValued(flag, giveIndex=True)
        if fv[0]:
            # flag is valued, use the index to return it
            return self.args[fv[1]]
        else:
            return None
            
    def isFlagValued(self, flag, giveIndex=False):
        # determine, if a flag is followed by a value other than another flag
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
    treeFiles = []
    def __init__(self, args):
        # sets the working directory of the script
        # if not set, then uses current directory
        self.args = args
        if self.args.isFlagValued('d'):
            # working directory is set by -d /path/to/wd
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
        # load tree files to work with from a file
        if self.fileExists(filename):
            with open(join(self.wd, filename), 'r') as file:
                # open the file to read
                self.treeFiles = []

                # process all lines
                for line in file.readlines():
                    # append and remove /r, /n
                    self.treeFiles += [line.rstrip()]

                countUnfiltered = len(self.treeFiles)
                # filter only existing files
                self.treeFiles = [f for f in self.treeFiles if self.fileExists(f)]
                
                deleted = countUnfiltered - len(self.treeFiles)
                if deleted > 0:
                    print(f"There were {deleted} non-existent files listed in '{filename}' and will be skipped")

        else:
            exit(f"Source file {filename} does not exist.")

    def loadFromDirectory(self, directory='.'):
        # load .fasta.tre file names form a defined directory relative to the working directory
        self.treeFiles = []

        for f in listdir(join(self.wd, directory)):
            # file is a file and has correct extension
            if isfile(join(self.wd, f)) and f[-10:] == TREE_FILE_EXTENSION:
                self.treeFiles += [f]

    def loadFilenames(self):
        if self.args.isFlagValued('f'):
            self.loadFromFile(self.args.getFlagValue('f'))
        else:
            self.loadFromDirectory('.')

class outputManager:
    def __init__(self, args, treeMgr):
        self.args = args
        if self.args.isFlagValued('o'):
            outfile = join(treeMgr.wd, self.args.getFlagValue('o'))
            if treeMgr.fileExists(outfile):
                if args.isFlagSet('r'):
                    # rewrite file
                    self.outfile = outfile
                    # delete contents
                    open(self.outfile, 'w').close()

                elif args.isFlagSet('c'):
                    # continuation will be attempted
                    self.outfile = outfile
                else:
                    # behavior unspecified, asking for flags
                    exit(f"Output file '{outfile}' already exists. Use -r to rewrite its contents or -c to continue progress")
            else:
                # file does not exist and will be created
                self.outfile = outfile
        else:
            # default output file
            self.outfile = join(treeMgr.wd, DEFAULT_OUTPUT_FILE)
                
class taxon:
    # node of a tree (taxon)
    def __init__(self, name, isRoot=False):
        self.name = name
        self.isRoot = isRoot
        # self.econdge = con
    # def setEdge(self, edge):
    #     if type(con) in [edge, node]:
    #         self.con = con
    #     else:
    #         exit(f"Failed attempt to set wrong type con {con} in {self}")

class node:
    def __init__(self):
        self.ways = [None, None, None]
        pass

    def getFork(self, edgeFrom):
        o = []
        for e in self.ways:
            if e != edgeFrom:
                o.append(e)
        if len(o) > 2:
            print(f"Node {self} called with wrong source node {edgeFrom}")
        else:
            return o

    def getLeft(self, edgeFrom):
        return self.getFork(edgeFrom)[0]

    def getRight(self, edgeFrom):
        return self.getFork(edgeFrom)[1]

class edge:
    a = None
    b = None
    bs = None
    # evaluated edge of a tree
    def __init__(self):
        pass

class phylogeneticTree:
    # parses fasta.tre file to build a tree, manages operations on it
    def __init__(self):
        self.rootElement = None
        self.edgeList = []
        self.taxonList = []
    
    def treFileLoad(self, treeMgr, path):
        # load tree file, strip header
        if not treeMgr.fileExists(path):
            exit(f"Tree file '{path}' could not be found")
        
        with open(path, 'r') as file:
            # load file
            # TODO remove headers
            tf = file.read().rstrip()   
                
class setDefinitions:
    # parses arguments to define subsets of the tree
    # one instance per tree
    def __init__(self, tree):
        pass
                
class criterionMatcher:
    # uses defined criteria to determine their validity on a subtree
    # one instance per defined criteria
    def __init__(self, tree, sets, criteria):
        pass

class definitionsManager:
    # parses all set and criteria definitions from arguments and builds their instances
    def __init__(self, args):
        pass




def main():
    # prepare argument parser
    args = argParser()
    # prepare input management according to flags
    treeMgr = treeFilesManager(args)
    # prepare output file according to flags
    outMgr = outputManager(args, treeMgr)
    print("Output will be:", outMgr.outfile)
    # load treeFiles according to flags
    treeMgr.loadFilenames()
    print(treeMgr.treeFiles)

    pt = phylogeneticTree()
    pt.treFileLoad(treeMgr, treeMgr.treeFiles[0])





if __name__ == "__main__":
    main()

def findMatchingBracket(data, i):
    status = 0
    o = 0
    for c in data[i+1:]:
        o += 1
        if status == 0 and c == ')':
            return i+o
        else:
            if c == ')':
                status -= 1
            elif c == '(':
                status += 1


taxons = []
edges = []

def recurse(data, con):
    newNode = node()

    newNode.ways[0] = con

    l = 0
    bs = None
    if data[l] == '(':
        r = findMatchingBracket(data, l)
        e = edge()
        edges.append(e)

        e.a = newNode
        e.b = recurse(data[l+1:r], e)
        addToNode = e
        bso = data[r+1:].find(',')
        bs = int(data[r+1:r+1+bso])
        e.bs = bs
        r += bso+2
        pass
    else:
        r = data[l:].find(',')
        addToNode = taxon(data[l:r])
        taxons.append(addToNode)
        r +=1

    x = data[r]

    newNode.ways[1] = addToNode

    l = r
    bs = None
    if data[l] == '(':
        r = findMatchingBracket(data, l)
        e = edge()
        edges.append(e)

        e.a = newNode
        e.b = recurse(data[l+1:r], e)
        addToNode = e
        bs = int(data[r+1:])
        e.bs = bs
        pass
    else:
        addToNode = taxon(data[l:])
        taxons.append(addToNode)


    newNode.ways[2] = addToNode

    return newNode



        







with open("phylo.fasta.tre", 'r') as tree:
    data = tree.read().rstrip()[1:-2]
    data = sub(r':[0-9]\.[0-9]*', '', data)


i = data.find(',')
root = taxon(data[:i], True)
taxons.append(root)

x = recurse(data[i+1:], root)

print(edges)
print(taxons)



pass


# recBrackets(data, 0, 'start')



