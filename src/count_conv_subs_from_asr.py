#!/usr/bin/env python3


"""Provides counts of convergent events per branch combo"""

import sys
import argparse
import tree_reader
from itertools import combinations
from collections import Counter


def read_subs(path):
    """takes as input a tab-delimited substitution file produced by
    summarise_asr_over_tree.py
    Consists of parent  descendant  subs[comma-separated, format A265V]"""
    subsDict = {}
    with open(path, "r") as inf:
        next(inf)  # skip header
        for line in inf.readlines():
            line = line.strip().split("\t")
            parentDescendant = (line[0], line[1])
            try:
                subsStr = line[2].split(",")
            except IndexError:  # no subs
                continue
            subs = [(x[0], int(x[1:-1]), x[-1]) for x in subsStr]
            subsDict[parentDescendant] = subs

    return subsDict


def get_ref_pos_dict(seqDict):
    allCorrDict = {}
    for k in seqDict.keys():
        corrDict = {}
        aln_pos = 1
        seq_pos = 1
        for i in seqDict[k]:
            if i == "-":
                corrDict[aln_pos] = 0
                aln_pos += 1
            else:
                corrDict[aln_pos] = seq_pos
                aln_pos += 1
                seq_pos += 1
        allCorrDict[k] = corrDict
    return allCorrDict


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="Tree file containing newick-formatted \
                        tree to run analysis on")
    parser.add_argument("subsfile", help="tsv with subs produced by \
                        summarise_asr_over_tree.py")
    parser.add_argument("branches", help="space-separated parent-daughter \
                        comparisons in format parent,daughter",
                        type=str, nargs="+")
    parser.add_argument("-a", "--atleast", help="sub must occur in at least \
                        n branches (default 2)", type=int, default=2)
    parser.add_argument("-nc", "--nonconv", help="also write \
                        non-site-overlapping changes on each branch",
                        action="store_true")
    args = parser.parse_args()

    branches = [(x.split(",")[0], x.split(",")[1]) for x in args.branches]

    subs_dict = read_subs(args.subsfile)

    curroot = next(tree_reader.read_tree_file_iter(args.tree))
    curroot.number_nodes()
    # make dict of tree nodes indexed by label
    node_dict = {}
    for n in curroot.iternodes():
        node_dict[n.label] = n
    
    branch_comb_dict = {}  # key is order, value is list of tuples
    i = 2
    while i <= args.atleast:
        branch_combs = combinations(branches, i)
        good_combs = []
        for b in branch_combs:
            sub_combs = combinations(b, 2)
            remove = False
            for s in sub_combs:
                # start arbitrarily with left parent
                n = node_dict[s[0][0]]
                going = True
                while going:
                    p = n.parent
                    if p.parent is None:  # reach root without encountering
                        going = False
                    if p.label == s[1][1]:  # right descendant is in parents
                        sys.stderr.write(f"skipping combination {b} as in "
                                         f"subcombination {s}, {s[1]} is "
                                         f"direct ancestor of {s[0]}\n")
                        remove = True
                        break
                    else:
                        n = p

                n = node_dict[s[0][1]]  # left descendant
                for c in n.iternodes():
                    if c.label == s[1][0]:  # right parent is in descendants
                        sys.stderr.write(f"skipping combination {b} as in "
                                         f"subcombination {s}, {s[1]} is "
                                         f"direct descendant of {s[0]}\n")
                        remove = True
                        break

            if not remove:  # no two-way combos are bad
                good_combs.append(b)
        branch_comb_dict[i] = good_combs
        i += 1

    # print(branch_comb_dict)

    out_dict = {}
    for order, combos in branch_comb_dict.items():
        for branch_comb in combos:  # specific branch combs at a given level
            out_dict[branch_comb] = {}
            subs = []
            for branch in branch_comb:
                try:
                    subs += [s for s in subs_dict[branch]]
                except KeyError:  # no subs on this branch
                    sys.stderr.write(f"no substitutions along {branch} "
                                     f"skipping {branch_comb}\n")
                    break
            pos = [s[1] for s in subs]
            if len(set(pos)) < len(subs):  # more than one branch has sub
                pos = [k for k, v in Counter(pos).items() if v == order]
                for p in pos:
                    end_states = [s[2] for s in subs if s[1] == p]
                    if len(set(end_states))  == 1:
                        try:
                            out_dict[branch_comb]["CONV"] += 1
                        except KeyError:
                            out_dict[branch_comb]["CONV"] = 1
                    else:
                        try:
                            out_dict[branch_comb]["DIV"] += 1
                        except KeyError:
                            out_dict[branch_comb]["DIV"] = 1
    
    print("order\tbranches\tconv\tdiv")
    for k, v in out_dict.items():
        order = len(k)
        try:
            conv = v["CONV"]
        except KeyError:
            conv = 0
        try:
            div = v["DIV"]
        except KeyError:
            div = 0
        branches = " ".join([",".join(i) for i in k])
        print(f"{order}\t{branches}\t{conv}\t{div}")

            

                


    # allSubsPos = []
    # for b in branchCombs:
    #     if args.groups is not None:
    #         contained = False
    #         for g in args.groups:
    #             g = [n for n in g.split(",")]
    #             if set(chain.from_iterable(b)).issubset(set(g)):
    #                 sys.stderr.write("skipping %s contained in %s\n"
    #                                     % (b, g))
    #                 contained = True
    #         if contained:
    #             continue
    #     subsPos = []
    #     for n in b:
    #         subsPos.append(set([x[1] for x in subs[n]]))
    #     subsInAll = set.intersection(*subsPos)
    #     allSubsPos.append(subsInAll)

    # allSubsPosUnion = sorted(list(set.union(*allSubsPos)))

    # for pos in allSubsPosUnion:
    #     endStates = []
    #     count_dict = {}
    #     for b in branches:
    #         endStates += [s[2] for s in subs[b] if pos == s[1]]
    #         if len(endStates) == 0:
    #             continue
    #     for b in branches:
    #         count_dict[b] = {}
    #         subsInPos = [s[0] + str(s[1]) + s[2] for s in subs[b]
    #                         if pos == s[1]]
    #         if len(subsInPos) == 0:
    #             continue
    #         # if args.groups is not None:
    #         #     uniq = len(branches)
    #         #     for g in args.groups:
    #         #         g = [n for n in g.split(",")]
    #         #         if set(chain.from_iterable(b)).issubset(set(g)):
    #         #             uniq -= 1
    #         #     if uniq < args.atleast:
    #         #         print("Uniq less than at least")
    #         #         continue
    #         endState = [s[2] for s in subs[b] if pos == s[1]][0]
    #         if endStates.count(endState) >= args.atleast:
    #             """make changes here in future to add requirements of 2-,
    #             3-way overlap etc."""
    #             try:
    #                 count_dict[b]["CONV"] += 1
    #             except KeyError:
    #                 count_dict[b]["CONV"] = 1
    #         else:
    #             try:
    #                 count_dict[b]["COIN"] += 1
    #             except KeyError:
    #                 count_dict[b]["COIN"] = 1
    # if args.nonconv:
    #     for b in branches:
    #         for s in subs[b]:
    #             pos = s[1]
    #             if pos not in allSubsPosUnion:
    #                 try:
    #                     count_dict[b]["NONC"] += 1
    #                 except KeyError:
    #                     count_dict[b]["NONC"] = 1
    # # subsPos = []
    # # for n in nodes:
    # #     subsPos.append(set([x[1] for x in subs[n]]))

    # # subsInAll = set.intersection(*subsPos)

    # # print(subsInAll)
    # print(count_dict)
