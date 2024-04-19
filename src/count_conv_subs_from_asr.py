#!/usr/bin/env python3


"""Provides counts of convergent events per branch combo"""

import sys
import argparse
from itertools import combinations
from calc_expected_conv import get_good_branch_combs, get_all_branches_labelled_tree
from summarise_asr_over_tree import parse_probs, add_subs, add_subs_robust
import sequence as sq
import newick as nwk


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
                subs_str = line[2].split(",")
            except IndexError:  # no subs
                continue
            substitutions = [(x[0], int(x[1:-1]), x[-1]) for x in subs_str]
            subs_dict[parent_descendant] = substitutions

    return subs_dict


def get_ref_pos_dict(seq_dict: dict) -> dict:
    """
    returns a dictionary showing the correspondence between each sequence position and alignment
    position for all sequences
    """
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


def count_conv_subs(comb: tuple[tuple], subs_dict: dict) -> dict:
    """
    function to count convergent and divergent substitutions in a branch combo of 
    ((par,desc ), (par, desc), ...) from a subs_dict produced by read_subs()
    """
    out_dict = {"CONV": 0, "DIV": 0}

    branch_subs = [s for br in comb for s in subs_dict[br]]
    # print(f"branch_subs: {branch_subs}")
    if all(branch_subs):
        positions = [sub[1] for sub in branch_subs]
        # print(f"positions: {positions}")
        duplicate_positions = set(p for p in positions if positions.count(p) > 1)

        for pos in duplicate_positions:
            end_states = {sub[2] for sub in branch_subs if sub[1] == pos}
            # print(f"end states: {end_states}")
            if len(end_states) == 1:
                out_dict["CONV"] += 1
            else:
                out_dict["DIV"] += 1
    else:
        raise ValueError(f"no substitutions along branches in {comb}")

    return out_dict


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="Tree file containing newick-formatted tree to run analysis \
                        on")
    parser.add_argument("alignment", help="FASTA-formatted alignment containing ancestral \
                        sequences with headers corresponding to node labels in tree")
    parser.add_argument("branches", help="space-separated parent-daughter comparisons in format \
                        parent,daughter", type=str, nargs="+")
    parser.add_argument("-a", "--all", help="get counts for all possible acceptable branch pairs, \
                        not just those in [branches ...]", action="store_true")
    parser.add_argument("-atl", "--atleast", help="sub must occur in at least n branches (default \
                        2)", type=int, default=2)
    parser.add_argument("-r", "--robust", help="count only substitutions where parent and child \
                        are unambiguous, i.e. PP > 0.8 (default False)", action="store_true")
    parser.add_argument("-p", "--probs", help="if robust. Tab-separated file of site, node, state,\
                        and probability, in the style of FastML's \
                        Ancestral_MaxMarginalProb_Char_Indel.txt")
    args = parser.parse_args()

    alignment = dict(sq.parse_fasta(args.alignment))

    curroot = nwk.parse_from_file(args.tree)
    curroot.number_nodes()
    # make dict of tree nodes indexed by label
    node_dict = {}
    for n in curroot.iternodes():
        node_dict[n.label] = n

    if args.all:
        branches_filt = []
        branches = get_all_branches_labelled_tree(curroot)
    else:
        branches = [(x.split(",")[0], x.split(",")[1]) for x in args.branches]

    subs = {br: [] for br in branches}

    branch_comb_dict = {}  # key is order, value is list of tuples
    o = 2
    while o <= args.atleast:
        branch_combs = combinations(branches, o)
        good_combs = get_good_branch_combs(branch_combs, curroot)
        branch_comb_dict[o] = good_combs
        o += 1

    if args.robust:
        if args.probs is None:
            sys.stderr.write("must specify a state probability file for robust mode\n")
            sys.exit()
        probs = parse_probs(args.probs)
        add_subs_robust(subs, alignment, probs)
    else:
        add_subs(subs, alignment)

    # print(branch_comb_dict)

    results = {}
    for order, combos in branch_comb_dict.items():
        for branch_comb in combos:  # specific branch combs at a given level
            try:
                results[branch_comb] = count_conv_subs(branch_comb, subs)
            except ValueError:
                sys.stderr.write(f"skipping combination {branch_comb} as no substitutions\n")
                continue

    print("order\tbranches\tconv\tdiv\tlength")
    for k, v in results.items():
        order = len(k)
        conv = v["CONV"]
        div = v["DIV"]
        try:
            length = v["path_length"]
        except KeyError:
            length = 0.
        branches = " ".join([",".join(i) for i in k])
        print(f"{order}\t{branches}\t{conv}\t{div}\t{length:.4f}")
