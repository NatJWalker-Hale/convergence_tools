#! /usr/bin/python3

import sys
import argparse
import tree_reader
import pandas as pd
from parse_fasta import parse_fasta
from treenode import Node


def parse_probs(file: str):
    df = pd.read_table(file, header=0)
    return df


def count_diffs(seq1: str, seq2: str, ignoreGap=True) -> list:
    """iterate over two aligned sequences and return differences
    as state1posstate2, e.g. A235V"""
    if len(seq1) != len(seq2):
        sys.stderr.write("sequences not equal in length!")
    rawDiffs = [(j[0], i+1, j[1]) for i, j in
                enumerate(zip(seq1, seq2)) if j[0] != j[1]]
    if ignoreGap:
        diffs = [j[0] + str(j[1]) + j[2] for j in rawDiffs if "-" not in j]
    else:
        diffs = [j[0] + str(j[1]) + j[2] for j in rawDiffs]
    return diffs


def get_anc_desc(node: Node) -> dict:
    """iterate over all nodes in a node-labelled tree and create
    a dictionary of ancestor descendant relationships, where each
    entry is one branch and value is empty list to be populated"""
    brDict = {}
    for n in node.iternodes(order="preorder"):
        if n.parent is None:
            continue  # skip root
        for c in n.children:
            brDict[(n.label, c.label, c.length)] = []
    return brDict


def add_subs(brDict, seqDict, ignoreGap=True):
    for k in brDict.keys():
        brDict[k] += count_diffs(seqDict[k[0]], seqDict[k[1]], ignoreGap)


def add_subs_robust(brDict, seqDict, probDF, ignoreGap=True):
    for k in brDict.keys():
        rawDiffs = count_diffs(seqDict[k[0]], seqDict[k[1]], ignoreGap)
        robDiffs = []
        for v in rawDiffs:
            p1 = probDF.loc[(probDF['Node'] == k[0]) & (probDF['Pos_on_MSA'] ==
                            int(v[1:-1]))]['CharProb'].iat[0]
            p2 = probDF.loc[(probDF['Node'] == k[1]) & (probDF['Pos_on_MSA'] ==
                            int(v[1:-1]))]['CharProb'].iat[0]
            if (p1 >= 0.8) and (p2 >= 0.8):
                robDiffs.append(v)
        brDict[k] += robDiffs


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="newick formatted tree, with ancestral \
                        node labels. Ignores root by default.")
    parser.add_argument("sequences", help="FASTA formatted \
                        alignment, including ancestral \
                        sequences")
    parser.add_argument("-g", "--gaps", help="include gaps (default False)",
                        type=bool, default=False)
    parser.add_argument("-r", "--robust", help="count only substitutions \
                        where parent and child are unambiguous, i.e. PP > 0.8 \
                        (default False)",
                        type=bool, default=False)
    parser.add_argument("-p", "--probs", help="if robust. Tab-separated file \
                        of site, node, state, and probability, in the style \
                        of FastML's Ancestral_MaxMarginalProb_Char_Indel.txt")
    args = parser.parse_args()

    if args.robust and args.probs is None:
        sys.stderr.write("must specify probability file (-p) to calculate " +
                         "robust substitutions\n")
        sys.exit()

    with open(args.tree, "r") as t:
        for s in t:
            s = s.strip()
            nwkString = s

    curroot = tree_reader.read_tree_string(nwkString)
    branches = get_anc_desc(curroot)
    #print(branches)

    seqs = dict([x for x in parse_fasta(args.sequences)])
    # print(seqs)

    if args.robust:
        probs = parse_probs(args.probs)
        add_subs_robust(branches, seqs, probs, not args.gaps)
    else:
        add_subs(branches, seqs, not args.gaps)

    print("parent\tchild\tsubs")
    # print(curroot.label)
    for k, v in branches.items():
        print(k[0] + "\t" + k[1] + "\t" + ",".join(v))
