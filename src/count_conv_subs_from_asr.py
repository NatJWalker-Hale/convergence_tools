#!/usr/bin/env python3


"""Provides counts of convergent events per branch combo"""

import sys
import argparse
from itertools import combinations
from collections import Counter
from calc_expected_conv import get_good_branch_combs
from calc_expected_conv import get_all_branches_labelled_tree
from calc_expected_conv import get_path_length_mrca
import tree_reader


def read_subs(path):
    """takes as input a tab-delimited substitution file produced by
    summarise_asr_over_tree.py
    Consists of parent  descendant  subs[comma-separated, format A265V]"""
    subs_dict = {}
    with open(path, "r", encoding='utf-8') as inf:
        next(inf)  # skip header
        for line in inf.readlines():
            line = line.strip().split("\t")
            parent_descendant = (line[0], line[1])
            try:
                subsStr = line[2].split(",")
            except IndexError:  # no subs
                continue
            subs = [(x[0], int(x[1:-1]), x[-1]) for x in subsStr]
            subs_dict[parent_descendant] = subs

    return subs_dict


def get_ref_pos_dict(seq_dict: dict) -> dict:
    all_corr_dict = {}
    for header in seq_dict.keys():
        corr_dict = {}
        aln_pos = 1
        seq_pos = 1
        for i in seq_dict[header]:
            if i == "-":
                corr_dict[aln_pos] = 0
                aln_pos += 1
            else:
                corr_dict[aln_pos] = seq_pos
                aln_pos += 1
                seq_pos += 1
        all_corr_dict[header] = corr_dict
    return all_corr_dict


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
    parser.add_argument("-a", "--all", help="Calculate expectation for all \
                        possible acceptable branch pairs, not just those in \
                        [branches ...]", action="store_true")
    parser.add_argument("-atl", "--atleast", help="sub must occur in at least \
                        n branches (default 2)", type=int, default=2)
    parser.add_argument("-nc", "--nonconv", help="also write \
                        non-site-overlapping changes on each branch",
                        action="store_true")
    args = parser.parse_args()

    subs = read_subs(args.subsfile)

    curroot = next(tree_reader.read_tree_file_iter(args.tree))
    curroot.number_nodes()
    # make dict of tree nodes indexed by label
    node_dict = {}
    for n in curroot.iternodes():
        node_dict[n.label] = n

    if args.all:
        branches_filt = []
        branches = get_all_branches_labelled_tree(curroot)
        for branch in branches:
            if branch in subs:
                branches_filt.append(branch)
        branches = branches_filt
    else:
        branches = [(x.split(",")[0], x.split(",")[1]) for x in args.branches]

    branch_comb_dict = {}  # key is order, value is list of tuples
    i = 2
    while i <= args.atleast:
        branch_combs = combinations(branches, i)
        good_combs = get_good_branch_combs(branch_combs, curroot)
        branch_comb_dict[i] = good_combs
        i += 1

    # print(branch_comb_dict)

    out_dict = {}
    for order, combos in branch_comb_dict.items():
        for branch_comb in combos:  # specific branch combs at a given level
            out_dict[branch_comb] = {}
            branch_subs = []
            for branch in branch_comb:
                try:
                    branch_subs += list(subs[branch])
                except KeyError:  # no subs on this branch
                    sys.stderr.write(f"no substitutions along {branch} "
                                     f"skipping {branch_comb}\n")
                    break
            pos = [s[1] for s in branch_subs]
            if len(set(pos)) < len(subs):  # more than one branch has sub
                pos = [k for k, v in Counter(pos).items() if v == order]
                for p in pos:
                    end_states = [s[2] for s in branch_subs if s[1] == p]
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
            if order == 2:
                lengths, _ = get_path_length_mrca(curroot, branch_comb[0][1],
                                                  branch_comb[1][1])
                out_dict[branch_comb]["path_length"] = sum(lengths)

    print("order\tbranches\tconv\tdiv\tlength")
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
        try:
            length = v["path_length"]
        except KeyError:
            length = 0
        branches = " ".join([",".join(i) for i in k])
        print(f"{order}\t{branches}\t{conv}\t{div}\t{length:.4f}")
