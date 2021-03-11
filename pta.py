#!/usr/bin/env python3

DEFAULT_OUTPUT_FILE = 'pta_output.csv' # for unspecified output file
TREE_FILE_EXTENSION = '.fasta.tre' # used when a whole directory of trees is set
RUN_QUIETLY = False # no prints will be executed
VERBOSE = False # enables some informational prints, error prints are unaffected

if RUN_QUIETLY:
    # this overrides the built-in print function to silence the whole program
    def overriden_print(*args):
        pass
    print = overriden_print

from sys import argv
from os import path, listdir
from os.path import isfile, join
import re

class argParser:
    def __init__(self, verbose=False):
        # load arguments from system interface
        # set basic variables within the instance
        self.args = argv[1:]
        self.fullArgsString = ' '.join(self.args)
        self.argsLength = len(self.args)

        if VERBOSE:
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
    def getOtherSide(self, e):
        if e in [self.a, self.b]:
            if e == self.a:
                return self.b
            else:
                return self.a
        else:
            exit("Wrong tree walk, programmer's fault.")

class phylogeneticTree:
    # parses fasta.tre file to build a tree, manages operations on it
    def __init__(self):
        self.root = None
        self.edges = []
        self.taxons = []
    
    def sortEdges(self):
        self.edges = sorted(self.edges, key=lambda e: -e.bs)

    def treFileLoad(self, treeMgr, path):
        # load tree file, strip header
        if not treeMgr.fileExists(path):
            exit(f"Tree file '{path}' could not be found")
        
        with open(path, 'r') as file:
            # load file
            # TODO remove headers
            data = file.read().rstrip()
            data = data.rstrip()[1:-2]
            data = re.sub(r':[0-9]\.[0-9]*', '', data)
            self.parseFile(data)

    def parseFile(self, data):
        i = data.find(',')
        self.root = taxon(data[:i], True)
        self.taxons.append(self.root)
        self.nodeRecurse(data[i+1:], self.root)

    def findMatchingBracket(self, data, i):
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


    def nodeRecurse(self, data, con):
        newNode = node()

        newNode.ways[0] = con

        l = 0
        bs = None
        if data[l] == '(':
            r = self.findMatchingBracket(data, l)
            bso = data[r+1:].find(',')
            bs = int(data[r+1:r+1+bso])

            e = edge()
            self.edges.append(e)

            e.a = newNode
            e.b = self.nodeRecurse(data[l+1:r], e)
            addToNode = e
            e.bs = bs
            r += bso+2
            pass
        else:
            r = data[l:].find(',')
            addToNode = taxon(data[l:r])
            self.taxons.append(addToNode)
            r +=1

        x = data[r]

        newNode.ways[1] = addToNode

        l = r
        bs = None
        if data[l] == '(':
            r = self.findMatchingBracket(data, l)
            e = edge()
            self.edges.append(e)

            e.a = newNode
            e.b = self.nodeRecurse(data[l+1:r], e)
            addToNode = e
            bs = int(data[r+1:])
            e.bs = bs
            pass
        else:
            addToNode = taxon(data[l:])
            self.taxons.append(addToNode)


        newNode.ways[2] = addToNode

        return newNode      
                
class setDefinitions:
    sets = {}
    # parses arguments to define subsets of the tree
    # one instance per tree
    def __init__(self, tree, args):
        self.argsString = args.getFlagValue('s')
        self.tree = tree

        sets = self.argsString.split(',')
        for name in sets:
            i = name.find('=')

            definition = re.sub(r'\*', '.*', name[i+1:] )
            if definition[:2] != '.*':
                definition = "^" + definition
            if definition[-2:] != '.*':
                definition = definition + "$"
            p = re.compile(definition)

            self.sets[name[:i]] = [ s.name for s in tree.taxons if p.match(s.name)]
        pass

    def countIntersects(self, taxons, setName):
        x = set(self.sets[setName])
        y = set(taxons)
        return length(x.intersection(y))

    def countDifference(self, taxons, setName):
        x = set(self.sets[setName])
        y = set(taxons)
        return length(x.difference(y))
               
class criteriaDefinitions:
    # uses defined criteria to determine their validity on a subtree
    # one instance per defined criteria
    crits = {}
    def __init__(self, tree, setDef, args):
        self.setDef = setDef
        self.tree = tree

        crits = args.getFlagValue('c').split(' ')
        for c in crits:
            i = c.find('=')
            self.crits[c[0:i]] = c[i+1:].split(',')
        if VERBOSE:
            print(self.crits) 
        pass

    def testCrit(self, c : dict, tl):
        if VERBOSE:
            print(f"Test criterium {c} = \"{self.crits[c]}\" on {tl}")

        for s in self.crits[c]:
            if s[-1:] in ['+','-']:
                maxCrit = True if s[-1] == '-' else False
                n = int(re.sub(r'^[A-Za-z_]', '', s[:-1]))
                sn = re.sub(r'[0-9]', '', s[:-1])

                diff = set(tl).intersection(self.setDef.sets[sn])
                diff.discard(self.tree.root.name)

                common = len(diff)
                if maxCrit:
                    if common > n:
                        return False
                else: # minCrit
                    if common < n:
                        return False
                
                pass

        unionCrit = set()
        for s in self.crits[c]:
            s = re.sub(r'[0-9]', '', s[:-1] if s[-1] in ['+','-'] else s)
                # normal sets
            unionCrit = unionCrit.union(self.setDef.sets[s])

        diff = set(tl).difference(unionCrit)
        diff.discard(self.tree.root.name)

        if len(diff) > 0:
            # found taxons from a sets complement
            return False



        return True

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
    if VERBOSE:
        print("Output will be:", outMgr.outfile)
    # load treeFiles according to flags
    treeMgr.loadFilenames()
    if VERBOSE:
        print(treeMgr.treeFiles)


    with open(outMgr.outfile, "w") as out:

        out.write(f"\"root_element\",")
        
        for c in args.getFlagValue('c').split(','):
            out.write(f"\"{c}\",")

        out.write(f"\"{args.fullArgsString}\"\n")

        for tr in treeMgr.treeFiles:

            pt = phylogeneticTree()
            pt.treFileLoad(treeMgr, tr)

            pt.sortEdges()

            # if VERBOSE:
            #     for i in pt.taxons:
            #         print("[R]" if i.isRoot else "   ", i.name)

            #     for i in pt.edges:
            #         print(f"Edge with BS: {i.bs} connects")
            #         for j in i.a.ways:
            #             if type(j) != edge:
            #                 print("\t\t", j.name)
            #             elif j != i:
            #                 print(f"\t\t another edge with BS: {j.bs}")
            #         print("\tand")
            #         for j in i.b.ways:
            #             if type(j) != edge:
            #                 print("\t\t", j.name)
            #             elif j != i:
            #                 print(f"\t\t another edge with BS: {j.bs}")
            #         print('')
            #     pass


            sd = setDefinitions(pt, args)

            cd = criteriaDefinitions(pt, sd, args)

            bestBoostraps = dict()

            out.write(f"\"{pt.root.name}\",")

            for ckey in cd.crits:
                bestBoostraps[ckey] = -1
                for e in pt.edges:
                    taxonList = []
                    subTrees(e.a, e, taxonList)
                    
                    if cd.testCrit(ckey, taxonList):
                        bestBoostraps[ckey] = e.bs
                        # print(f"Best bootstrap for {ckey} = {cd.crits[ckey]} is {e.bs}")
                        out.write(f"\"{e.bs}\",")
                        break
                if bestBoostraps[ckey] < 0:
                    # print(f"No valid bootstrap for {ckey} = {cd.crits[ckey]} found")
                    out.write(f"\"0\",")
                out.write("\"0\"\n")
                    


        
def subTrees(n, e, tl):
    for i in n.getFork(e):
        if type(i) == taxon:
            tl.append(i.name)
        elif type(i) == edge:
            subTrees(i.getOtherSide(n), i, tl)








   






if __name__ == "__main__":
    main()








# recBrackets(data, 0, 'start')



