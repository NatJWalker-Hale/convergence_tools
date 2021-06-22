#! /usr/bin/python3

import sys
import re
import argparse
import tree_reader


def repl(m):
    n = m.group(0).split(".")[0]
    return n


def label_for_diffsel(root, scenarioDict, nCond, outPath="diffsel_tree.nwk"):
    for n in root.iternodes(order="preorder"):
        if nCond == 2:
            if n.number in scenarioDict[1]:
                n.length = 1.0
            else:
                n.length = 0.0
        elif nCond == 3:
            if n.number in scenarioDict[1]:
                n.length = 1.0
            elif n.number in scenarioDict[2]:
                n.length = 2.0
            else:
                n.length = 0.0

    newNwkString = root.get_newick_repr(showbl=True)
    # now we need to strip the branch lengths to integers
    p = re.compile(":[0-9].[0-9]*")

    newNwkString = p.sub(repl, newNwkString)

    with open(outPath, "w") as outTree:
        outTree.write(newNwkString+";\n")


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("treeFile", help="tree file in newick format")
    parser.add_argument("scenario", help="text file containing convergent \
                        scenario in PCOC format. Can have up to two lines \
                        defining the foreground branches for 2 foreground \
                        conditions.")
    args = parser.parse_args()

    with open(args.treeFile, "r") as inf:
        nwkString = inf.readline().strip()
        curroot = tree_reader.read_tree_string(nwkString)

    curroot.number_tree()

    scenarios = {}
    with open(args.scenario, "r") as inf:
        nCond = 1
        for s in inf:
            if s != "\n":
                scenarios[nCond] = []
                for i in s.strip().split("/"):
                    scenarios[nCond] += [int(x) for x in i.split(",")]
                nCond += 1

    label_for_diffsel(curroot, scenarios, nCond)
