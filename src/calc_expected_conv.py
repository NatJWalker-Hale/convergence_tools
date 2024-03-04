#!/usr/bin/env python3


"""Calculate expected convergence and divergence given a tree and some
reconstructed states following Zou & Zhang (2015)"""


import sys
import argparse
from collections import Counter
from itertools import combinations, product
import numpy as np
import tree_reader
from treenode import Node
from parse_fasta import parse_fasta
from discrete_models import Discrete_model
from optimise_raxmlng import get_optimised_freqs


AAIDX = {"A": 0, "R": 1, "N": 2, "D": 3, "C": 4, "Q": 5, "E": 6, "G": 7,
         "H": 8, "I": 9, "L": 10, "K": 11, "M": 12, "F": 13, "P": 14, "S": 15,
         "T": 16, "W": 17, "Y": 18, "V": 19}
AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F',
      'P', 'S', 'T', 'W', 'Y', 'V']


def get_seq_dict_from_tree(tree: Node, seq_dict: dict) -> dict:
    """takes a tree and sequence dictionary and returns a dictionary of
    sequences that are tips in the tree"""
    out = {}
    for n in tree.iternodes():
        if n.istip:
            try:
                out[n.label] = seq_dict[n.label]
            except KeyError:
                sys.stderr.write(f"{n.label} not in provided FASTA - "
                                    f"extant sequences need to be "
                                    f"included\n")
                sys.exit()
    return out


def get_columns(seq_dict: dict) -> dict:
    """takes a dictionary of an alignment (key: name, value: sequence),
    and returns a dictionary of columns 
    (key: position, value: dict{name: state})"""
    col_dict = {}
    pos = 0
    for header, seq in seq_dict.items():
        for char in seq:
            try:
                col_dict[pos][header] = char
            except KeyError:
                col_dict[pos] = {}
                col_dict[pos][header] = char
            pos += 1
        pos = 0
    return col_dict


def get_site_specific_frequencies(col_dict: dict) -> dict:
    """takes a dictionary from get_columns and returns a dictionary of 
    site-specific frequencies"""
    freq_dict = {}  # key is pos, value is np array of freqs
    for pos, column in col_dict.items():
        freqs = []
        tot = 0
        counts = Counter([c for c in column.values()])
        tot += sum(x[1] for x in counts.items() if x[0] in AA)
        for state in AA:
            freqs.append(counts[state] / tot)
        freqs = np.array(freqs)
        freq_dict[pos] = freqs
    return freq_dict


def get_opt_site_specific_frequencies(col_dict: dict, tree: str) -> dict:
    """takes a dictionary from get_columns, writes individual column FASTAs,
    optimises ML frequencies with raxml-ng and returns the frequencies as a 
    dictionary"""
    out = {}
    for pos, col in col_dict.items():
        with open(f"tmp_{pos}.aln", "w", encoding='utf-8') as sitef:
            for header, char in col.items():
                sitef.write(f">{header}\n")
                sitef.write(f"{char}\n")
        freqs = get_optimised_freqs(tree, f"tmp_{pos}.aln", "JTT+FO")
        out[pos] = freqs
    with open("site_frequencies.tsv", "w", encoding="utf-8") as siteff:
        for pos, freq in out.items():
            siteff.write(f"{pos}\t{','.join([str(f) for f in freq])}\n")
    return out


def parse_site_frequencies_file(inf: str) -> dict:
    """parse a TSV of pos\tfreqs, where freqs is comma-separated"""
    freqs = {}
    with open(inf, "r", encoding="utf-8") as sff:
        for line in sff:
            line = line.strip().split("\t")
            freqs[int(line[0])] = np.array([float(x) for
                                            x in line[1].split(",")])
    return freqs


def parse_paml_rates(path: str) -> dict:
    """takes path to a paml 'rates' file containing estimated posterior mean 
    rate per site and returns a dictionary of {site: rate}"""
    out = {}
    with open(path, "r", encoding='utf-8') as inf:
        going = False
        reading = False
        for line in inf:
            if line.strip().startswith("Site"):
                going = True
                continue
            if line == "\n" and going and not reading:
                reading = True
                continue
            if line == "\n" and going and reading:
                break

            if going and reading:
                line = line.strip().split()
                out[int(line[0])] = float(line[-2])
    return out


def get_node_dict_from_tree(tree: Node) -> dict:
    """takes a tree with node labels and returns a dictionary of 
    {label: Node}"""
    out = {}
    for n in tree.iternodes():
        if not n.istip and n.label is None:
            sys.stderr.write("encountered internal node with no label. Might "
                             "want to run tree.number_nodes() first?\n")
            sys.exit()
        else:
            out[n.label] = n
    return out


def get_good_branch_combs(combs: list[tuple[tuple]],
                          tree: Node) -> list[tuple[tuple]]:
    """takes an input list of tuples of two tuples of length two (e.g.
    [(("N7", "N8"), ("N41", "N42")), (("N171", "N224"), ("N8", "N9'))]) 
    and a tree with node labels matching labels in tuples and determines 
    which combinations to keep"""
    out = []
    node_dict = get_node_dict_from_tree(tree)
    for c in combs:
        remove = False
        par1, ch1, par2, ch2 = (x for i in c for x in i)
        if par1 == par2:  # sisters
            sys.stderr.write(f"skipping combination {c} as {c[0]} is "
                             f"sister to {c[1]}\n")
            remove = True
        # start arbitrarily with left parent
        n = node_dict[par1]
        going = True
        while going:
            p = n.parent
            if p.parent is None:  # reach root without encountering
                going = False
            if p.label == ch2:  # right descendant is in parents
                sys.stderr.write(f"skipping combination {c} as {c[1]} is "
                                 f"direct ancestor of {c[0]}\n")
                remove = True
            else:
                n = p

        n = node_dict[ch1]  # left descendant
        for n1 in n.iternodes():
            if n1.label == par2:  # right parent is in descendants
                sys.stderr.write(f"skipping combination {c} as {c[1]} is "
                                 f"direct descendant of {c[0]}\n")
                remove = True

        if not remove:
            out.append(c)
    return out


def get_mrca(tree: Node, node1: str, node2: str) -> Node:
    """returns the node that is the mrca for two labelled nodes"""
    for n in tree.iternodes():
        if n.label == node1:
            going = True
            while going:
                p = n
                if node2 in [i.label for i in p.iternodes()]:
                    going = False
                else:
                    n = p.parent
    return p


def map_state_array(aa_char: str) -> np.array:
    """returns a dim20 array of zeros with idx aa = 1"""
    idx = AAIDX[aa_char]
    array = np.zeros(20)
    array[idx] = 1.
    return array


def get_path_length_mrca(tree: Node, node1: str, node2: str) -> tuple:
    """returns the sum of the branch lengths between two labelled nodes and
    their mrca"""
    mrca = None
    lengths = []
    for n in tree.iternodes():
        if n.label == node1:
            going = True
            brlen = 0.
            while going:
                if node2 in [i.label for i in n.iternodes()]:  # node IS mrca
                    going = False
                    mrca = n
                else:
                    brlen += n.length
                    n = n.parent
    lengths.append(brlen)
    for n in tree.iternodes():
        if n.label == node2:
            going = True
            brlen = 0.
            while going:
                if node1 in [i.label for i in n.iternodes()]:  # node IS mrca
                    going = False
                    mrca = n
                else:
                    brlen += n.length
                    n = n.parent
    lengths.append(brlen)
    return lengths, mrca


def get_path_length_root(tree: Node, node1: str):
    """returns the sum of branch lengths between the node and the root"""
    for n in tree.iternodes():
        if n.label == node1:
            going = True
            brlen = 0.
            while going:
                p = n
                brlen += n.length
                if p.parent is None:  # root
                    going = False
                else:
                    n = p.parent
    return brlen


def main(combs: list[tuple[tuple]], tree: Node, model: Discrete_model,
         ancestor_columns: dict, rates: dict = None,
         site_frequencies: dict = None) -> dict:
    """calculates expected convergence and divergence for a given tree and
    model, conditioning on ancestral states in ancestor_columns"""
    out = {}
    node_dict = get_node_dict_from_tree(tree)
    for c in combs:
        par1, ch1, par2, ch2 = (x for i in c for x in i)
        conv_prob_sum = 0.
        div_prob_sum = 0.
        lengths, mrca = get_path_length_mrca(tree, par1, par2)
        if mrca.label in (par1, par2):
            # if either of the nodes is the MRCA, we calculate probs for
            # this node by conditioning on the state one node back, to
            # avoid fixed states
            lengths = [x + mrca.length for x in lengths]
            mrca = mrca.parent
        # get lengths between node and MRCA
        sites = 0
        if rates is None:  # do these once and save time
            p0 = model.get_P(lengths[0])
            p1 = model.get_P(lengths[1])
            p01 = model.get_P(node_dict[ch1].length)
            p11 = model.get_P(node_dict[ch2].length)
        for pos, col in ancestor_columns.items():   # going by site
            map_state = col[mrca.label]  # get MAP state of MRCA
            if map_state == "-":
                continue  # skip cases where mrca MAP is gap
            nodes_in = [x for i in c for x in i]
            states = [col[n] for n in nodes_in]
            if "-" in states:
                continue  # skip cases where any MAPs are inferred as gaps
            if site_frequencies is not None:
                # set up model for site specific freqs
                model = Discrete_model()
                model.set_rate_JTT()
                model.set_frequencies(site_frequencies[pos])
                with np.errstate(divide="raise"):
                    # account for cases where 0 freqs cause problems
                    try:
                        model.scale_rate_matrix()
                    except FloatingPointError:
                        sys.stderr.write(f"skipping site {pos} where "
                                         f"sparse frequencies caused a "
                                         f"problem\n")
                        continue
            if rates is not None:
                p0 = model.get_P(lengths[0], rates[pos])
                p1 = model.get_P(lengths[1], rates[pos])
                p01 = model.get_P(node_dict[ch1].length,
                                  rates[pos])
                p11 = model.get_P(node_dict[ch2].length,
                                  rates[pos])
            else:
                p0 = model.get_P(lengths[0])
                p1 = model.get_P(lengths[1])
                p01 = model.get_P(node_dict[ch1].length)
                p11 = model.get_P(node_dict[ch2].length)
            # dict for probs
            sites += 1
            node_probs = {}
            # probs for parents
            map_array = map_state_array(map_state)
            node_probs[par1] = np.matmul(map_array, p0)  #1x20
            node_probs[par2] = np.matmul(map_array, p1)  #1x20
            # conditional probs for children
            conditionals_l = np.zeros((20,20))
            conditionals_r = np.zeros((20,20))
            for i in range(20):
                map_array = np.zeros(20)
                map_array[i] = 1.
                prob_l = np.matmul(map_array, p01)  # 1x20
                conditionals_l[i] = prob_l
                prob_r = np.matmul(map_array, p11)  # 1x20
                conditionals_r[i] = prob_r
                # conditionals now holds the probability of desc being
                # in state j (col) given starting in state i (row)
            joint_probs_l = np.matmul(np.diag(node_probs[par1]),
                                        conditionals_l)
            joint_probs_r = np.matmul(np.diag(node_probs[par2]),
                                        conditionals_r)
            for perm in product(range(20), repeat=4):
                i, j, k, l = perm
                if i != j and k != l and j == l:
                    # print(f"Branch 1: {AA[i]} -> {AA[j]}")
                    # print(f"Branch 2: {AA[k]} -> {AA[l]}")
                    conv_prob_sum += (joint_probs_l[i, j] *
                                      joint_probs_r[k, l])
                elif i != j and l != k and j != l:
                    div_prob_sum += (joint_probs_l[i, j] *
                                     joint_probs_r[k, l])
        out[c] = [sites, conv_prob_sum, div_prob_sum]
    return out


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="Tree file containing newick-formatted \
                        tree to run analysis on")
    parser.add_argument("ancestors", help="FASTA-formatted alignment of MAP \
                        sequences for nodes, with names matching labels \
                        produced by --label")
    parser.add_argument("branches", help="space-separated parent-daughter \
                        comparisons in format parent,daughter (at least two)",
                        type=str, nargs="+")
    group = parser.add_mutually_exclusive_group(required = False)
    group.add_argument("-f", "--empirical_frequencies", help="Use \
                        frequencies calculated from extant sequences. \
                        'ancestors' should contain sequences for both \
                        nodes and tips", action="store_true")
    group.add_argument("-sf", "--site_frequencies", help="Use per-site \
                        frequencies calculated from extant sequences. \
                        'ancestors' should contain sequences for both \
                        nodes and tips", action="store_true")
    parser.add_argument("-sff", "--site_frequencies_file", help="if using \
                        site frequencies, read from file. Format: tab \
                        separated, first column position (0-indexed) \
                        second column comma separated frequencies (20)")
    parser.add_argument("-r", "--rates", help="File containing estimated \
                        posterior mean rate for each site from PAML")
    parser.add_argument("-d", "--divergent", help="Also calculate the \
                        probability of divergent substitutions",
                        action="store_true")
    parser.add_argument("-l", "--label", help="Print labelled tree and \
                        quit", action="store_true")
    args = parser.parse_args()
    if args.site_frequencies_file:
        args.site_frequencies = True

    curroot = next(tree_reader.read_tree_file_iter(args.tree))

    curroot.number_nodes()
    # make dict of tree nodes indexed by label
    nodes = get_node_dict_from_tree(curroot)

    if args.label:
        print(f"{curroot.get_newick_repr(False, False)};")
        sys.exit()

    if args.rates:
        rates_dict = parse_paml_rates(args.rates)
    else:
        rates_dict = None

    ancs = dict(parse_fasta(args.ancestors))
    anc_cols = get_columns(ancs)

    # set up model
    modJTT = Discrete_model()
    modJTT.set_rate_JTT()
    if args.empirical_frequencies:
        extants = get_seq_dict_from_tree(curroot, ancs)
        modJTT.calc_empirical_freqs(extants)
        modJTT.set_empirical_freqs()
    else:
        modJTT.set_model_freqs()
    modJTT.scale_rate_matrix()

    if args.site_frequencies:
        if args.site_frequencies_file:
            site_freqs = parse_site_frequencies_file(args.site_frequencies_file)
        else:
            extants = get_seq_dict_from_tree(curroot, ancs)
            extant_cols = get_columns(extants)
            site_freqs = get_opt_site_specific_frequencies(extant_cols,
                                                           args.tree)
    else:
        site_freqs = None

    branches = [(x.split(",")[0], x.split(",")[1]) for x in args.branches]
    if len(branches) == 1:
        print("More than one branch required")
        sys.exit()

    if len(branches) == 2:
        branch_combs = [tuple(branches)]
    else:
        branch_combs = list(combinations(branches, 2))

    good_combs = get_good_branch_combs(branch_combs, curroot)

    results = main(combs=good_combs, tree=curroot, model=modJTT,
                   ancestor_columns=anc_cols, rates=rates_dict,
                   site_frequencies=site_freqs)

    print("branch_comb\tsites\tconv\tdiv")
    for comb, res in results.items():
        print(f"{comb}\t{res[0]}\t{res[1]}\t{res[2]}")
