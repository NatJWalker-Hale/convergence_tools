#! /usr/bin/python3

import re
import sys
import argparse
import tree_reader
from parse_fasta import parse_fasta


def get_group(label):
    if "_" in label:
        p = re.compile("^.*?_")
        group = p.search(label).group(0)[:-1]
        return group
    else:
        return label


def label_for_tdg09(root, seqDict, scenarioDict, cond0, cond1,
                    outTreePath="tdg09_tree.nwk",
                    outAlnPath="tdg09_aln.phy"):
    newSeqDict = {}
    for n in root.iternodes(order="preorder"):
        if n.number in scenarioDict[1]:
            if n.istip:
                n.label = n.label.replace("_" + str(n.number),
                                          "", 1)
                newLab = cond1 + "_" + n.label
                newSeqDict[newLab] = seqDict[n.label]
                n.label = newLab
            else:
                n.label = cond1
        else:
            if n.istip:
                n.label = n.label.replace("_" + str(n.number),
                                          "", 1)
                newLab = cond0 + "_" + n.label
                newSeqDict[newLab] = seqDict[n.label]
                n.label = newLab
            else:
                n.label = cond0

    for n in root.iternodes(order="preorder"):
        # print(n.label)
        for c in n.children:
            # print(c.label)
            nGroup = get_group(n.label)
            cGroup = get_group(c.label)
            if nGroup != cGroup:
                n.label = re.sub(nGroup, nGroup + "_GS", n.label)

    with open(outTreePath, "w") as outTree:
        outTree.write(root.get_newick_repr(showbl=True)+";\n")

    with open(outAlnPath, "w") as outAln:
        outAln.write(str(len(newSeqDict.keys())) + " " +
                     str(len(list(newSeqDict.values())[0])) + "\n")
        for k, v in newSeqDict.items():
            outAln.write(k + "\t" + v + "\n")


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
    curroot.number_tree()
    seqs = dict([x for x in parse_fasta(args.alignment)])

    scenarios = {}
    with open(args.scenario, "r") as inf:
        nCond = 1
        for s in inf:
            if s != "\n":
                scenarios[nCond] = []
                for i in s.strip().split("/"):
                    scenarios[nCond] += [int(x) for x in i.split(",")]
                nCond += 1

    # print(scenario)

    label_for_tdg09(curroot, seqs, scenarios, args.cond0, args.cond1)
