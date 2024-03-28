#!/usr/bin/env python3


"""
script to infer substitutions from reconstructed ancestral sequences over a tree
"""


import sys
import argparse
import pandas as pd
from parse_fasta import parse_fasta
from phylo import Node
import newick as nwk


def parse_probs(file: str):
    """
    read state posterior probabilities from a CSV-formatted file
    """
    df = pd.read_table(file, header=0)
    return df


def count_diffs(seq1: str, seq2: str, gaps=False) -> list:
    """iterate over two aligned sequences and return differences
    as state1posstate2, e.g. A235V"""
    if len(seq1) != len(seq2):
        raise ValueError("sequences not equal in length!")
    raw_diffs = [(j[0], i+1, j[1]) for i, j in
                enumerate(zip(seq1, seq2)) if j[0] != j[1]]
    if not gaps:
        diffs = [j for j in raw_diffs if "-" not in j]
    else:
        diffs = list(raw_diffs)
    return diffs


def get_anc_desc(node: Node) -> dict:
    """iterate over all nodes in a node-labelled tree and create a dictionary of ancestor descendant 
    relationships, where each entry is one branch and value is empty list to be populated"""
    br_dict = {}
    for n in node.iternodes():
        for c in n.children:
            br_dict[(n.label, c.label)] = []
    return br_dict


def add_subs(br_dict: dict, seq_dict: dict, gaps: bool=False):
    """
    populates the individual elements of a list of branches with substitution information drawn
    from reconstructed and extant sequences in seq_dict
    """
    for br in br_dict:
        try:
            br_dict[br] = count_diffs(seq_dict[br[0]], seq_dict[br[1]], gaps)
        except KeyError as e:
            raise KeyError("sequence not found") from e


def add_subs_robust(br_dict: dict, seq_dict: dict, prob_df: pd.array, gaps: bool=False,
                    thresh: float=0.8):
    """
    populates the individual elements of a list of branches with substitution information drawn
    from reconstructed and extant sequences in seq_dict, first determining from a file of state
    probabilities if both parent and descendant states have > 
    """
    for br in br_dict:
        try:
            raw_diffs = count_diffs(seq_dict[br[0]], seq_dict[br[1]], gaps)
            rob_diffs = []
            for d in raw_diffs:
                p1 = prob_df.loc[(prob_df['Node'] == br[0]) & (prob_df['Pos_on_MSA'] ==
                                 d[1])]['CharProb'].iat[0]
                p2 = prob_df.loc[(prob_df['Node'] == br[1]) & (prob_df['Pos_on_MSA'] ==
                                 d[1])]['CharProb'].iat[0]
                if (p1 >= thresh) and (p2 >= thresh):
                    rob_diffs.append(d)
            br_dict[br] = rob_diffs
        except KeyError as e:
            raise KeyError("sequence not found") from e


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="newick formatted tree, with ancestral \
                        node labels (will label if not present). Ignores root by default.")
    parser.add_argument("sequences", help="FASTA formatted \
                        alignment, including ancestral \
                        sequences")
    parser.add_argument("-g", "--gaps", help="include gaps",
                        action="store_true")
    parser.add_argument("-r", "--robust", help="count only substitutions \
                        where parent and child are unambiguous, i.e. PP > 0.8 \
                        (default False)",
                        type=bool, default=False)
    parser.add_argument("-p", "--probs", help="if robust. Tab-separated file of site, node, state,\
                        and probability, in the style of FastML's \
                        Ancestral_MaxMarginalProb_Char_Indel.txt")
    args = parser.parse_args()

    if args.robust and args.probs is None:
        sys.stderr.write("must specify probability file (-p) to calculate " +
                         "robust substitutions\n")
        sys.exit()

    curroot = nwk.parse_from_file(args.tree)
    if "" in [n.label for n in curroot.iternodes() if not n.istip]:
        curroot.number_nodes()
        sys.stderr.write("Here is your tree with labelled nodes. Ancestral sequence headers should "
                         "match these labels\n")
        sys.stderr.write(f"{curroot.to_string};\n")
    branches = get_anc_desc(curroot)
    #print(branches)

    seqs = dict(parse_fasta(args.sequences))
    # print(seqs)
    incl_gaps = bool(args.gaps)

    if args.robust:
        probs = parse_probs(args.probs)
        add_subs_robust(branches, seqs, probs, incl_gaps)
    else:
        add_subs(branches, seqs, incl_gaps)

    # print(branches)
    print("parent\tchild\tsubs")
    # print(curroot.label)
    for branch, subs in branches.items():
        print(f"{branch[0]}\t{branch[1]}\t{','.join([''.join(map(str, x)) for x in subs])}")
    