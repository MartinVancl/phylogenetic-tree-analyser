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
        # create full arguments list to prefix output file
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
                # sets working directory
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
                # counts how many listed filenames are non-existent
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
            # get filenames from external file
            self.loadFromFile(self.args.getFlagValue('f'))
        else:
            # use current directory as default
            self.loadFromDirectory('.')

class outputManager:
    # manages exporting to csv
    def __init__(self, args, treeMgr):
        self.args = args
        if self.args.isFlagValued('o'):
            # output file is set and is not default
            outfile = join(treeMgr.wd, self.args.getFlagValue('o'))
            if treeMgr.fileExists(outfile):
                # checks existence of the output file
                if args.isFlagSet('r'):
                    # rewrite file
                    self.outfile = outfile
                    # delete contents
                    open(self.outfile, 'w').close()

                elif args.isFlagSet('a'):
                    # continuation will be attempted
                    self.outfile = outfile
                else:
                    # behavior unspecified, asking for flags
                    exit(f"Output file '{outfile}' already exists. Use -r to rewrite its contents or -a to continue progress")
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
    # connects three edges or taxons
    def __init__(self):
        self.ways = [None, None, None]
        pass

    def getFork(self, edgeFrom):
        # gives you the other two connected elements, than one passed as edgeFrom
        o = []
        for e in self.ways:
            if e != edgeFrom:
                o.append(e)
        if len(o) > 2:
            # errors if tree walking asked fork from non-present edge/taxon
            print(f"Node {self} called with wrong source node {edgeFrom}")
        else:
            return o

    def getLeft(self, edgeFrom):
        # gives first of forks
        return self.getFork(edgeFrom)[0]

    def getRight(self, edgeFrom):
        # gives second of forks
        return self.getFork(edgeFrom)[1]

class edge:
    # evaluated edge of a tree
    def __init__(self):
        self.a = None
        self.b = None
        self.bs = None
        pass
    def getOtherSide(self, e):
        # returns the other side of the edge
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
        # order edgelist by its bootstrap value to trim calculations
        self.edges = sorted(self.edges, key=lambda e: -e.bs)

    def treFileLoad(self, treeMgr, path):
        # load tree file, strip header
        if not treeMgr.fileExists(path):
            exit(f"Tree file '{path}' could not be found")
        
        with open(path, 'r') as file:
            # load file
            data = file.read().rstrip()
            # remove encapsulation bracket pair
            data = data.rstrip()[1:-2]
            # remove unimportant edge evaluations (length)
            data = re.sub(r':[0-9]\.[0-9]*', '', data)
            self.parseFile(data)

    def parseFile(self, data):
        # start building tree from data
        # get root taxon name from the start
        i = data.find(',')
        self.root = taxon(data[:i], True)
        self.taxons.append(self.root)
        # run recursive tree builder
        self.nodeRecurse(data[i+1:], self.root)

    def findMatchingBracket(self, data, i):
        # state automat for finding corresponding right bracket position
        status = 0
        o = 0
        for c in data[i+1:]:
            # go through string by characters
            o += 1
            if status == 0 and c == ')':
                # matching states and bracket found
                return i+o
            else:
                if c == ')':
                    # reduce state
                    status -= 1
                elif c == '(':
                    # increase state
                    status += 1


    def nodeRecurse(self, data, con):
        # create a new node to build around
        newNode = node()
        # link "parent" element to this node
        newNode.ways[0] = con

        # go through data string and find corresponding two elements
        l = 0
        bs = None
        if data[l] == '(':
            # starts with subtree, find it's end bracked and recursively process
            r = self.findMatchingBracket(data, l)
            bso = data[r+1:].find(',')
            bs = int(data[r+1:r+1+bso])
            # create connecting edge
            e = edge()
            self.edges.append(e)
            # connect the edge to current node
            e.a = newNode
            # connect the edge to new subtree
            e.b = self.nodeRecurse(data[l+1:r], e)
            addToNode = e
            e.bs = bs
            # increase position in data to process next element
            r += bso+2
            pass
        else:
            # starts with taxon (leaf), create and append the taxon
            r = data[l:].find(',')
            addToNode = taxon(data[l:r])
            self.taxons.append(addToNode)
            r +=1

        # append new link to this node
        newNode.ways[1] = addToNode

        # continue to the other node of this fork
        l = r
        bs = None
        if data[l] == '(':
            # element is a subtree, process it
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
            # element is a taxon, create and connect it
            addToNode = taxon(data[l:])
            self.taxons.append(addToNode)


        newNode.ways[2] = addToNode
        # new node was created and will be returned to be connected to parent node
        return newNode      
                
class setDefinitions:
    sets = {}
    # parses arguments to define subsets of the tree
    # one instance per tree
    def __init__(self, tree, args):
        self.argsString = args.getFlagValue('s')
        self.tree = tree
        # parse all argument-defined sets
        sets = self.argsString.split(',')
        for name in sets:
            # get definition and name split index
            i = name.find('=')
            # get definition string and replace wildcards with regex equivalent
            definition = re.sub(r'\*', '.*', name[i+1:] )
            if definition[:2] != '.*':
                # add string start regex
                definition = "^" + definition
            if definition[-2:] != '.*':
                # add string end regex
                definition = definition + "$"
            # build the definition's regex
            p = re.compile(definition)
            # create this definition's list of valid taxons
            self.sets[name[:i]] = [ s.name for s in tree.taxons if p.match(s.name)]
        pass

    def countIntersects(self, taxons, setName):
        # number of intersections with all taxons
        x = set(self.sets[setName])
        y = set(taxons)
        return length(x.intersection(y))

    def countDifference(self, taxons, setName):
        # number of differences with all taxons
        x = set(self.sets[setName])
        y = set(taxons)
        return length(x.difference(y))
               
class criteriaDefinitions:
    # uses defined criteria to determine their validity on a subtree
    # one instance per defined criteria
    crits = {}
    def __init__(self, setDef, args):
        self.tree = setDef.tree
        self.setDef = setDef
        # parse criteria from arguments
        crits = args.getFlagValue('c').split(' ')
        for c in crits:
            # split name and definition
            i = c.find('=')
            # add all set citerie of this name as list
            self.crits[c[0:i]] = c[i+1:].split(',')
        if VERBOSE:
            print(self.crits) 
        pass

    def testCrit(self, c : dict, tl):
        if VERBOSE:
            print(f"Test criterium {c} = \"{self.crits[c]}\" on {tl}")

        for s in self.crits[c]:
            if s[-1:] in ['+','-']:
                # if criterium is has number-constraing
                # set which constraint it is (max or min)
                maxCrit = True if s[-1] == '-' else False
                # parse set name for the constraint
                n = int(re.sub(r'^[A-Za-z_]', '', s[:-1]))
                # parse the number for the constraint
                sn = re.sub(r'[0-9]', '', s[:-1])
                # get difference of all tested taxons and this set
                diff = set(tl).intersection(self.setDef.sets[sn])
                # remove root taxon
                diff.discard(self.tree.root.name)
                # count how many taxons are left
                common = len(diff)
                if maxCrit:
                    # max criterim test
                    if common > n:
                        # failed test
                        return False
                else:
                    # min criterium test
                    if common < n:
                        # failed test
                        return False
                
                pass

        unionCrit = set()
        for s in self.crits[c]:
            # testing "clean content" criteria
            # removes number and +/- from definitions where they were
            s = re.sub(r'[0-9]', '', s[:-1] if s[-1] in ['+','-'] else s)
            # now we have "pretend" normal sets
            # append (union) the set to a big set of all
            unionCrit = unionCrit.union(self.setDef.sets[s])
        # differ tested taxons and the built set
        diff = set(tl).difference(unionCrit)
        # remove root taxons (if it was part of defined set, it would be removed already)
        diff.discard(self.tree.root.name)
        # calculate how many taxons are left
        if len(diff) > 0:
            # found taxons from a sets complement, so they are not just from it
            return False

        return True

class definitionsManager:
    # parses all set and criteria definitions from arguments and builds their instances
    def __init__(self, args):
        pass

def subTrees(n, e, tl):
    for i in n.getFork(e):
        if type(i) == taxon:
            tl.append(i.name)
        elif type(i) == edge:
            subTrees(i.getOtherSide(n), i, tl)

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

    # prepare file for export
    with open(outMgr.outfile, "w") as out:
        # start writing header
        out.write(f"\"root_element\",")
        # write criteria definitions header
        for c in args.getFlagValue('c').split(','):
            out.write(f"\"{c}\",")
        # add whole arguments as part of a header
        out.write(f"\"{args.fullArgsString}\"\n")


        # run processing for all tree files
        for tr in treeMgr.treeFiles:
            # create new phylogenetic tree for current file
            pt = phylogeneticTree()
            pt.treFileLoad(treeMgr, tr)

            # sort edges
            # sort edges to trim tests
            pt.sortEdges()

            # build set definitions
            sd = setDefinitions(pt, args)
            # build criteria definitions
            cd = criteriaDefinitions(sd, args)
            
            # starting from highest bootstrap
            # gives you best result as soon as criteria pass
            if VERBOSE:
                # print taxons and edges and what they connect
                for i in pt.taxons:
                    # print taxons
                    print("[R]" if i.isRoot else "   ", i.name)

                for i in pt.edges:
                    # print edges
                    print(f"Edge with BS: {i.bs} connects")
                    for j in i.a.ways:
                        # print what it connects from
                        if type(j) != edge:
                            print("\t\t", j.name)
                        elif j != i:
                            print(f"\t\t another edge with BS: {j.bs}")
                    print("\tand")
                    for j in i.b.ways:
                        # print what it connects to
                        if type(j) != edge:
                            print("\t\t", j.name)
                        elif j != i:
                            print(f"\t\t another edge with BS: {j.bs}")
                    print('')
                pass

            
            bestBoostraps = dict()

            # start writing row, write out name 
            out.write(f"\"{pt.root.name}\",")

            # run criteria test for each one
            for ckey in cd.crits:
                # set impossibly bad bootstrap
                bestBoostraps[ckey] = -1
                for e in pt.edges:
                    # go through all sorted edges
                    taxonList = []
                    # build subtree in direction to root split by this edge
                    subTrees(e.a, e, taxonList)
                    
                    if cd.testCrit(ckey, taxonList):
                        # criterium passed, set it as best
                        bestBoostraps[ckey] = e.bs
                        # best criterium found and will be exported
                        out.write(f"\"{e.bs}\",")
                        break
                if bestBoostraps[ckey] < 0:
                    # no valid subtree found, exporting as bootstrap 0
                    out.write(f"\"0\",")
                # append ending thiny so i don't have to care about commas
                out.write("\"0\"\n")
                    


        








   






if __name__ == "__main__":
    main()








# recBrackets(data, 0, 'start')



