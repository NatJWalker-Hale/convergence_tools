#! /usr/bin/python3

import re
import sys
import argparse
import tree_reader
from copy import deepcopy
from parse_fasta import parse_fasta


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="newick tree with PAML-style labels to \
                        transform to tdg09 format")
    parser.add_argument("alignment", help="FASTA-formatted alignment matching \
                        sequences in tree")
    parser.add_argument("scenario", help="PCOC-formatted scenario in txt file")
    parser.add_argument("cond0", help="label for background condition")
    parser.add_argument("cond1", help="label for foreground condition")
    args = parser.parse_args()

    curroot = [x for x in tree_reader.read_tree_file_iter(args.tree)][0]
    seqDict = dict([x for x in parse_fasta(args.alignment)])
    newSeqDict = {}
    curroot.number_tree()

    scenario = []
    with open(args.scenario, "r") as inf:
        for line in inf:
            for i in line.strip().split("/"):
                scenario += [int(x) for x in i.split(",")]
    # print(scenario)

    for n in curroot.iternodes(order="preorder"):
        if n.number in scenario:
            if n.istip:
                n.label = n.label.replace("_" + str(n.number),
                                          "", 1)
                n.label += "_" + args.cond1
            else:
                n.label = args.cond1
        else:
            if n.istip:
                n.label = n.label.replace("_" + str(n.number),
                                          "", 1)
                n.label += "_" + args.cond0
            else:
                n.label = args.cond0
    for n in curroot.iternodes(order="preorder"):
        for c in n.children:
            if n.label[:-2] != c.label[:-2]:
                n.label += "_" + "GS"
             
    # for n in curroot.iternodes(order="preorder"):
    #     if n.label == "#1":  # only two conditions allowed
    #         n.label = ""  # strip label
    #         for leaf in n.leaves_fancy():
    #             newLab = args.cond1 + "_" + leaf.label
    #             newSeqDict[newLab] = seqDict[leaf.label]
    #             leaf.label = newLab
    #     if n.istip:
    #         if "#1" in n.label:
    #             newLab = args.cond1 + "_" + n.label[:-2]
    #             newSeqDict[newLab] = seqDict[n.label[:-2]]
    #             n.label = newLab
    #         elif n.label not in newSeqDict.keys():
    #             newLab = args.cond0 + "_" + n.label
    #             newSeqDict[newLab] = seqDict[n.label]
    #             n.label = newLab

    with open("tdg09_tree.nwk", "w") as outTree:
        outTree.write(curroot.get_newick_repr(True)+";\n")

    with open("tdg09_aln.phy", "w") as outAln:
        outAln.write(str(len(newSeqDict.keys())) + " " +
                     str(len(list(newSeqDict.values())[0])) + "\n")
        for k, v in newSeqDict.items():
            outAln.write(k + "\t" + v + "\n")
  
