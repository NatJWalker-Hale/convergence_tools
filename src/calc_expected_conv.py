#!/usr/bin/env python3


"""
Calculate expected convergence and divergence given a tree and some reconstructed states following 
Zou & Zhang (2015)
"""


import sys
import argparse
from itertools import combinations
import numpy as np
import sequence as sq
from newick import parse_from_file
from phylo import Node, getMRCATraverse
from discrete_models import DiscreteModel
from optimise_raxmlng import get_optimised_freqs


def get_seq_dict_from_tree(tree: Node, seq_dict: dict) -> dict:
    """takes a tree and sequence dictionary and returns a dictionary of
    sequences that are tips in the tree"""
    out = {}
    for n in tree.iternodes():
        if n.istip:
            try:
                out[n.label] = seq_dict[n.label]
            except KeyError as e:
                raise KeyError(f"{n.label} not in provided FASTA - extant "
                               f"sequences need to be included") from e
    return out


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
            siteff.write(f"{pos}\t" + ','.join([f"{num:.6f}" for num in freq]) + "\n")
    return out


def get_node_dict_from_tree(tree: Node) -> dict:
    """
    takes a tree with node labels and returns a dictionary of {label: Node}
    """
    out = {}
    for n in tree.iternodes():
        if not n.istip and n.label is None:
            raise ValueError("encountered internal node with no label. Might "
                             "want to run tree.number_nodes() first?\n")
        out[n.label] = n
    return out


def get_all_branches_labelled_tree(tree: Node) -> list[tuple]:
    """
    take a node-labelled tree and returns all branches as a list of 2-tuples of (parent-label, 
    descendant-label), not including root or branches immediately descending root
    """
    out = []
    for n in tree.iternodes(order="postorder"):
        par = n.parent
        if par is None:  # skip root
            continue
        if par.parent is None:  # skip children of root
            continue
        if None in (n.label, par.label):
            raise ValueError("found node with no labels - this only works on "
                             "labelled trees. Run tree.number_nodes() and "
                             "retry")
        out.append((par.label, n.label))
    return out

# refactored with Claude's assistance
def get_good_branch_combs(combs: list[tuple[tuple]],
                          tree: Node) -> list[tuple[tuple]]:
    """takes an input list of tuples of two tuples of length two (e.g.
    [(("N7", "N8"), ("N41", "N42")), (("N171", "N224"), ("N8", "N9'))]) 
    and a tree with node labels matching labels in tuples and determines 
    which combinations to keep"""
    out = []
    node_dict = get_node_dict_from_tree(tree)
    ancestors = {node: set(node.get_ancestors(True)) for node in
                 node_dict.values()}
    for c in combs:
        par1, ch1, par2, ch2 = (node_dict[x] for x in (x for i in c for x in i))

        if par1 == par2:  # sisters
            sys.stderr.write(f"skipping combination {c} as {c[0]} is "
                             f"sister to {c[1]}\n")
            continue
        # of course, if branch 1 is an ancestor of branch 2 then branch 2 is a
        # descendant of branch 1 and vice versa
        if ch1 in ancestors[par2]:  # direct ancestor
            sys.stderr.write(f"skipping combination {c} as {c[0]} is "
                             f"direct ancestor of {c[1]}\n")
            continue

        if ch2 in ancestors[par1]:
            sys.stderr.write(f"skipping combination {c} as {c[0]} is "
                             f"direct descendant of {c[1]}\n")
            continue

        out.append(c)

    return out


# deprecate
# def get_mrca(tree: Node, node1: str, node2: str) -> Node:
#     """returns the node that is the mrca for two labelled nodes"""
#     for n in tree.iternodes():
#         if n.label == node1:
#             going = True
#             while going:
#                 p = n
#                 if node2 in [i.label for i in p.iternodes()]:
#                     going = False
#                 else:
#                     n = p.parent
#     return p


def map_state_array(aa_char: str) -> np.array:
    """returns a dim20 array of zeros with idx aa = 1"""
    aaidx = {"A": 0, "R": 1, "N": 2, "D": 3, "C": 4, "Q": 5, "E": 6, "G": 7,
         "H": 8, "I": 9, "L": 10, "K": 11, "M": 12, "F": 13, "P": 14, "S": 15,
         "T": 16, "W": 17, "Y": 18, "V": 19}
    idx = aaidx[aa_char]
    array = np.zeros(20)
    array[idx] = 1.
    return array


def get_path_length_mrca(node1: Node, node2: Node) -> tuple:
    """
    for two nodes, returns a list of summed lengths between themselves and the
    MRCA, not including the MRCA branch length
    """
    lengths = []
    mrca = getMRCATraverse(node1, node2)
    lengths.append(sum([l.length for l in
                        node1.get_ancestors(True, mrca)][:-1]))
    lengths.append(sum([l.length for l in
                        node2.get_ancestors(True, mrca)][:-1]))
    return lengths, mrca


def get_path_length_root(node: Node):
    """returns the sum of branch lengths between the node and the root"""
    return sum(a.length for a in node.get_ancestors(True))


def compute_transition_probs(model: DiscreteModel, brlens: list[float],
                             rate: float=1.0):
    """
    computes the four transition probability matrices necessary for joint
    convergence probability calculation for two branches
    """
    p0 = model.get_P(brlens[0], rate)
    p1 = model.get_P(brlens[1], rate)
    p01 = model.get_P(brlens[2], rate)
    p11 = model.get_P(brlens[3], rate)
    return p0, p1, p01, p11


def compute_transition_probs_site_freqs(frequencies: np.array,
                                        brlens: list[float], rate: float=1.0):
    """
    computes the four transition probability matrices necessary for joint
    convergence probability calculation for two branches, rescaling the model
    matrix to use site-specific frequencies
    """
    model = DiscreteModel()
    model.set_rate_JTT()
    model.set_frequencies(frequencies)
    with np.errstate(divide="raise"):
        # account for cases where 0 freqs cause problems
        try:
            model.scale_rate_matrix()
        except FloatingPointError as exc:
            raise FloatingPointError("sparse frequencies") from exc
    p0 = model.get_P(brlens[0], rate)
    p1 = model.get_P(brlens[1], rate)
    p01 = model.get_P(brlens[2], rate)
    p11 = model.get_P(brlens[3], rate)
    return p0, p1, p01, p11


def main(combs: list[tuple[tuple]], tree: Node, model: DiscreteModel,
         ancestor_columns: dict, rates: dict = None,
         site_frequencies: dict = None, div = False) -> dict:
    """
    calculates expected convergence and divergence for a given tree and
    model, conditioning on ancestral states in ancestor_columns
    """
    out = {}
    node_dict = get_node_dict_from_tree(tree)
    for c in combs:
        par1, ch1, par2, ch2 = (x for i in c for x in i)
        conv_prob_sum = 0.
        div_prob_sum = 0.
        lengths, mrca = get_path_length_mrca(node_dict[par1], node_dict[par2])
        if mrca.label in (par1, par2):
            # if either of the nodes is the MRCA, we calculate probs for
            # this node by conditioning on the state one node back, to
            # avoid fixed states
            lengths = [x + mrca.length for x in lengths]
            mrca = mrca.parent
        # get lengths between node and MRCA
        lengths += [node_dict[x].length for x in (ch1, ch2)]
        # add child lengths to lengths
        sites = 0
        if rates is None:  # do these once and save time
            p0, p1, p01, p11 = compute_transition_probs(model, lengths)
        for pos, col in ancestor_columns.items():   # going by site
            rate = 1.0
            map_state = col[mrca.label]  # get MAP state of MRCA
            if map_state == "-":
                continue  # skip cases where mrca MAP is gap
            nodes_in = [x for i in c for x in i]
            states = [col[n] for n in nodes_in]
            if "-" in states:
                continue  # skip cases where any MAPs are inferred as gaps
            if site_frequencies is not None:
                # set up model for site specific freqs
                if rates is not None:
                    rate = rates[pos]
                p0, p1, p01, p11 = compute_transition_probs_site_freqs(
                    frequencies=site_frequencies[pos], brlens=lengths,
                    rate=rate
                    )
            else:
                if rates is not None:
                    rate = rates[pos]
                    p0, p1, p01, p11 = compute_transition_probs(model,
                                                                brlens=lengths,
                                                                rate=rate)
            # dict for probs
            sites += 1
            # probs for parents
            map_array = map_state_array(map_state)
            par1_probs = np.diag(map_array @ p0)  #1x20
            par2_probs = np.diag(map_array @ p1)  #1x20
            # conditional probs for children
            # conditionals_l = np.zeros((20,20))
            # conditionals_r = np.zeros((20,20))
            idmat = np.eye(20)
            conditionals_l = idmat @ p01
            conditionals_r = idmat @ p11
            # for i in range(20):
            #     map_array = np.zeros(20)
            #     map_array[i] = 1.
            #     prob_l = np.matmul(map_array, p01)  # 1x20
            #     conditionals_l[i] = prob_l
            #     prob_r = np.matmul(map_array, p11)  # 1x20
            #     conditionals_r[i] = prob_r
            #     # conditionals now holds the probability of desc being
            #     # in state j (col) given starting in state i (row)
            joint_probs_l = par1_probs @ conditionals_l
            joint_probs_r = par2_probs @ conditionals_r
            # for perm in product(range(20), repeat=4):
            #     i, j, k, l = perm
            #     if i != j and k != l and j == l:
            #         # print(f"Branch 1: {AA[i]} -> {AA[j]}")
            #         # print(f"Branch 2: {AA[k]} -> {AA[l]}")
            #         conv_prob_sum += (joint_probs_l[i, j] *
            #                           joint_probs_r[k, l])
            #     if div:
            #         if i != j and l != k and j != l:
            #             div_prob_sum += (joint_probs_l[i, j] *
            #                             joint_probs_r[k, l])
            # help from claude to figure out meshgrid
            ids = np.array(np.meshgrid(*[np.arange(20)] * 4, indexing='ij')).reshape(4, -1).T
            # ids.shape: (160000, 4)
            cond1 = np.logical_and.reduce([ids[:, 0] != ids[:, 1],
                                           ids[:, 2] != ids[:, 3],
                                           ids[:, 1] == ids[:, 3]])
            # cond1.shape: (160000,)
            # np.nonzero(cond1)[0].shape: (7220,)
            # that is, the 7220 cases satisfying the constraint (i != j and k != l and i==l)
            # here we are getting the logical vectors corresponding to conditions specified
            # in the three-way conditional statement. cond1 then corresponds to which entries
            # in ids satisfy the requirement
            conv_prob_sum += np.sum(joint_probs_l[ids[cond1, 0], ids[cond1, 1]] *
                                    joint_probs_r[ids[cond1, 2], ids[cond1, 3]])
            # by indexing by ids[cond1, 0], we are getting the rows of joint_probs_l that
            # match the constraint for each branch's states, vectorised
            # so that in the first contribution, we are indexing
            # joint_probs_l[0, 1] * joint_probs_r[0, 1], which is the joint probability of
            # convergent substitutions from AA[0] (A) to AA[1] (R)
            if div:
                cond2 = np.logical_and.reduce([ids[:, 0] != ids[:, 1],
                                               ids[:, 2] != ids[:, 3],
                                               ids[:, 1] != ids[:, 3]])     
                div_prob_sum += np.sum(joint_probs_l[ids[cond2, 0], ids[cond2, 1]] *
                                       joint_probs_r[ids[cond2, 2], ids[cond2, 3]])
            # this numpy approach provides _dramatic_ speedups
        child_lengths, _ = get_path_length_mrca(node_dict[ch1],
                                                node_dict[ch2])
        sys.stderr.write(f"{c}\t{sites}\t{conv_prob_sum:.4f}\t"
                         f"{div_prob_sum:.4f}\t{sum(child_lengths):.4f}\n")
        out[c] = [sites, conv_prob_sum, div_prob_sum, sum(child_lengths)]
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
    parser.add_argument("-a", "--all", help="Calculate expectation for all \
                        possible acceptable branch pairs, not just those in \
                        [branches ...]", action="store_true")
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

    curroot = parse_from_file(args.tree)

    curroot.number_nodes()
    # make dict of tree nodes indexed by label
    nodes = get_node_dict_from_tree(curroot)

    if args.label:
        print(f"{curroot.to_string()};")
        sys.exit()

    # print(get_path_length_mrca(curroot, "N224", "N172"))
    # print(get_path_length_mrca_new(nodes["N224"], nodes["N172"]))
    # print([l.label for l in nodes["N224"].get_ancestors(True, nodes["N171"])])
    # sys.exit()

    if args.rates:
        rates_dict = sq.parse_paml_rates(args.rates)
    else:
        rates_dict = None

    ancs = dict(sq.parse_fasta(args.ancestors))
    anc_cols = sq.get_columns(ancs)

    # set up model
    modJTT = DiscreteModel()
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
            site_freqs = sq.parse_site_frequencies_file(args.site_frequencies_file)
        else:
            extants = get_seq_dict_from_tree(curroot, ancs)
            site_freqs = sq.get_site_specific_frequencies(extants)
            with open("site_frequencies.tsv", "w", encoding="utf-8") as sff:
                for pos, freq in site_freqs.items():
                    sff.write(f"{pos}\t")
                    sff.write(",".join([f"{f:.4f}" for f in freq]) + "\n")
    else:
        site_freqs = None

    if args.all:
        branches = get_all_branches_labelled_tree(curroot)
    else:
        branches = [(x.split(",")[0], x.split(",")[1]) for x in args.branches]
    if len(branches) == 1:
        sys.stderr.write("More than one branch required\n")
        sys.exit(1)

    if len(branches) == 2:
        branch_combs = [tuple(branches)]
    else:
        branch_combs = list(combinations(branches, 2))

    good_combs = get_good_branch_combs(branch_combs, curroot)

    results = main(combs=good_combs, tree=curroot, model=modJTT,
                   ancestor_columns=anc_cols, rates=rates_dict,
                   site_frequencies=site_freqs, div=args.divergent)

    print("branch_comb\tsites\tconv\tdiv\tlength")
    for comb, res in results.items():
        branches = " ".join([",".join(b) for b in comb])
        print(f"{branches}\t{res[0]}\t{res[1]:.4f}\t{res[2]:.4f}\t"
              f"{res[3]:.4f}")
