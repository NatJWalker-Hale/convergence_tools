#! /usr/bin/python3

import sys
import argparse
import tree_reader


def prep_for_msd(root, scenario):
    for n in curroot.iternodes(order="preorder"):
        if n.istip:
            if n.number in scenario:
                stateDict[n.label] = 1
            else:
                stateDict[n.label] = 0

    with open("msd_tree.nwk", "w") as outTree:
        outTree.write(curroot.get_newick_repr(showbl=True)+";\n")

    with open("msd_states.tsv", "w") as outFile:
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

    scenario = []
    with open(args.scenario, "r") as inf:
        for line in inf:
            for i in line.strip().split("/"):
                scenario += [int(x) for x in i.split(",")]

    prep_for_msd(curroot, scenario)
