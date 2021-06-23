#! /usr/bin/python3

import sys
import argparse
import tree_reader


def prep_for_msd(root, scenarioDict, outTreePath="msd_tree.nwk",
                 outStatePath="msd_states.tsv"):
    stateDict = {}

    for n in root.iternodes(order="preorder"):
        if n.istip:
            if n.number in scenarioDict[1]:
                stateDict[n.label] = 1
            else:
                stateDict[n.label] = 0

    with open(outTreePath, "w") as outTree:
        outTree.write(root.get_newick_repr(showbl=True)+";\n")

    with open(outStatePath, "w") as outFile:
        for k, v in stateDict.items():
            outFile.write(k + "\t" + str(v) + "\n")


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="newick tree")
    parser.add_argument("scenario", help="PCOC-formatted scenario in txt file")
    args = parser.parse_args()

    curroot = [x for x in tree_reader.read_tree_file_iter(args.tree)][0]
    stateDict = {}

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

    prep_for_msd(curroot, scenario)
