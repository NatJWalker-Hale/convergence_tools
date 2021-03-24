#! /usr/bin/python3

import sys
import argparse
import tree_reader


def parse_cond_file(path):
    condDict = {}
    with open(path, "r") as inf:
        for line in inf:
            line = line.strip().split("\t")
            condDict[line[0]] = line[1]
    return condDict

if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("conditionFile", help="tab separated file of \
                        taxon\tcondition, one per line")
    parser.add_argument("treeFile", help="tree file in newick format with \
                        matching taxa")
    args = parser.parse_args()

    with open(args.treeFile, "r") as inf:
        nwkString = inf.readline().strip()
        curroot = tree_reader.read_tree_string(nwkString)

    cond = parse_cond_file(args.conditionFile)

    for n in curroot.iternodes(order="preorder"):
        descStates = [cond[x] for x in n.lvsnms()]
        if all(descStates):  # node is mono for cond
            n.label = int(set(descStates)[0])
        else:
            n.label = 0

    print(curroot)
