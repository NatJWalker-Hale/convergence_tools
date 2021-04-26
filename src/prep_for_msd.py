#! /usr/bin/python3

import sys
import argparse
import tree_reader


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="newick tree with PAML-style labels to \
                        get msd states")
    args = parser.parse_args()

    curroot = [x for x in tree_reader.read_tree_file_iter(args.tree)][0]
    stateDict = {}

    for n in curroot.iternodes(order="preorder"):
        if n.label == "#1":  # only two conditions allowed
            n.label = ""  # strip label
            for leaf in n.leaves_fancy():
                stateDict[leaf.label] = 1
        if n.istip:
            if "#1" in n.label:
                n.label = n.label[:-2]
                stateDict[n.label] = 1
            elif n.label not in stateDict.keys():
                stateDict[n.label] = 0

    with open("msd_tree.nwk", "w") as outTree:
        outTree.write(curroot.get_newick_repr(True)+";\n")

    with open("msd_states.tsv", "w") as outFile:
        for k, v in stateDict.items():
            outFile.write(k + "\t" + str(v) + "\n")
