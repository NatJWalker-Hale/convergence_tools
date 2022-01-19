#! /usr/bin/python3

import sys
import argparse
import tree_reader
from collections import Counter
from scipy.stats import chi2
from scipy.stats import multinomial
from parse_fasta import parse_fasta


def get_columns(seqDict):
    colDict = {}
    pos = 0
    for k, v in seqDict.items():
        for i in v:
            try:
                colDict[pos].append(i)
            except KeyError:
                colDict[pos] = []
                colDict[pos].append(i)
            pos += 1
        pos = 0
    return colDict


def get_col_counts(colDict):
    colCountDict = {}
    aa = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F",
          "P", "S", "T", "W", "Y", "V"]
    for k, v in colDict.items():
        counts = []
        for i in aa:
            try:
                c = Counter(v)[i]
            except KeyError:
                c = 0
            counts.append(c)
        colCountDict[k] = counts
    return colCountDict


def smooth_prop_from_counts(countDict):
    """following simplexFromCounts in Rey et al.'s calc_multinomial.
    reapportions 1% of total density to each element with a count of
    zero."""
    smoothPropDict = {}
    for k, v in countDict.items():
        v2 = [x * 100 for x in v]
        for i in range(len(v2)):
            if v2[i] == 0:
                v2[i] = 1
        tot = sum(v2)
        for i in range(len(v2)):
            v2[i] = v2[i] / tot
        smoothPropDict[k] = (v, v2)
    return smoothPropDict


def calc_multinomial(smoothPropDict1, smoothPropDict2, smoothPropDictAll):
    """following multinomial_lrt from Rey et al.'s calc_multinomial"""
    resDict = {}
    for k in smoothPropDict1.keys():
        rv1 = multinomial(sum(smoothPropDict1[k][0]), smoothPropDict1[k][1])
        rv2 = multinomial(sum(smoothPropDict2[k][0]), smoothPropDict2[k][1])
        rvAll = multinomial(sum(smoothPropDictAll[k][0]),
                            smoothPropDictAll[k][1])
        print(k)
        lk1 = rv1.logpmf(smoothPropDict1[k][0])
        print(lk1)
        lk2 = rv2.logpmf(smoothPropDict2[k][0])
        print(lk2)
        lkAll = rvAll.logpmf(smoothPropDictAll[k][0])
        print(lkAll)
        lr = 2 * (-lkAll - lk1 + lk2)
        lrt = chi2.sf(lr, 1)
        resDict[k] = (lr, lrt)
    return resDict


def calc_col_prop(colDict):
    colPropDict = {}
    aa = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F",
          "P", "S", "T", "W", "Y", "V"]
    for k, v in colDict.items():
        tot = sum([x[1] for x in Counter(v).items() if x[0] != "-"])
        props = []
        for i in aa:
            try:
                c = Counter(v)[i]
            except KeyError:
                c = 0
            props.append(c / tot)
            colPropDict[k] = props
    return colPropDict


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("aln", help="alignment, in fasta")
    parser.add_argument("tree", help="tree, in newick")
    parser.add_argument("scenario", help="text file containing scenario(s), \
                        in PCOC format")
    args = parser.parse_args()

    curroot = [x for x in tree_reader.read_tree_file_iter(args.tree)][0]
    curroot.number_tree()
    seqs = dict([x for x in parse_fasta(args.aln)])

    scenarios = {}
    with open(args.scenario, "r") as inf:
        nCond = 1
        for s in inf:
            if s != "\n":
                scenarios[nCond] = []
                for i in s.strip().split("/"):
                    scenarios[nCond] += [int(x) for x in i.split(",")]
                nCond += 1
    scenario = scenarios[1]

    foreground = {}
    background = {}

    for n in curroot.iternodes(order="preorder"):
        if n.istip:
            if n.number in scenario:
                foreground[n.label] = seqs[n.label]
            else:
                background[n.label] = seqs[n.label]

    foreCols = get_columns(foreground)
    backCols = get_columns(background)
    allCols = get_columns(seqs)

    foreCounts = get_col_counts(foreCols)
    backCounts = get_col_counts(backCols)
    allCounts = get_col_counts(allCols)

    foreProps = smooth_prop_from_counts(foreCounts)
    backProps = smooth_prop_from_counts(backCounts)
    allProps = smooth_prop_from_counts(allCounts)
    # print(allProps)

    results = calc_multinomial(foreProps, backProps, allProps)

    print("pos\tlr\tlrt")
    for k, v in results.items():
        print(str(k) + "\t" + str(v[0]) + "\t" + str(v[1]))
