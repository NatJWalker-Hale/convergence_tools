#!/usr/bin/env python3


"""
bootstrap expected number of convergences and divergences for a given alignment and tree using 
simulation
"""


import sys
import argparse
from itertools import combinations
import sequence as sq
import newick as nwk
from phylo import Node
from sim_site_specific_freqs import sim_alignment_gene_freqs, sim_alignment_site_freqs
from sim_site_specific_freqs import sim_sites_indelible
from calc_expected_conv import get_good_branch_combs, get_node_dict_from_tree
from count_conv_subs_from_asr import count_conv_subs
from summarise_asr_over_tree import add_subs


def run_sim_gf(tree: Node, alignment: dict, branch_combs: list[tuple[tuple]], model: str="JTT",
               gamma: tuple[int, float]=None, use_anc: int=1, reps: int=100,
               write: bool=False) -> dict:
    """
    function to run and parse expect convergence and divergence by simulation, using gene-wide
    frequencies. Frequencies can be specified or if not will be calculated from the input alignment
    """
    results_dict = {comb: {"ex_conv": 0, "ex_div": 0, "obs_conv": 0, "obs_div": 0, "sims_c_lt": 0,
                           "sims_c_gt": 0, "sims_d_lt": 0, "sims_d_gt": 0, "p_c": 0., "p_d": 0.,
                           "r_c": 0., "r_d": 0., "reps": 0} for comb in branch_combs}
    ex_br_dict = {br: [] for x in branch_combs for br in x}
    obs_br_dict = {br: [] for x in branch_combs for br in x}
    add_subs(obs_br_dict, alignment)
    for comb in branch_combs:
        try:
            obs_res = count_conv_subs(comb, obs_br_dict)
            results_dict[comb]["obs_conv"] = obs_res["CONV"]
            results_dict[comb]["obs_div"] = obs_res["DIV"]
        except ValueError:
            sys.stderr.write(f"skipping combination {comb} as no subs on at least one branch\n")
            continue
    i = 0
    while i < reps:
        sim_aln = sim_alignment_gene_freqs(tree, alignment, model, gamma, use_anc)
        sim_aln = sq.insert_gaps_by_seq(alignment, sim_aln)
        if write:
            sq.write_fasta(sim_aln, f"rep_{i}.pep.fa")
        add_subs(ex_br_dict, sim_aln)

        for comb in branch_combs:
            try:
                ex_res = count_conv_subs(comb, ex_br_dict)
                results_dict[comb]["ex_conv"] += ex_res["CONV"]
                results_dict[comb]["ex_div"] += ex_res["DIV"]
                if ex_res["CONV"] >= results_dict[comb]["obs_conv"]:
                    results_dict[comb]["sims_c_gt"] += 1
                if ex_res["CONV"] <= results_dict[comb]["obs_conv"]:
                    results_dict[comb]["sims_c_lt"] += 1
                if ex_res["DIV"] >= results_dict[comb]["obs_div"]:
                    results_dict[comb]["sims_d_gt"] += 1
                if ex_res["DIV"] <= results_dict[comb]["obs_div"]:
                    results_dict[comb]["sims_d_lt"] += 1
                results_dict[comb]["reps"] += 1
            except ValueError:
                sys.stderr.write(f"skipping combination {comb} as no subs on at least on branch\n")
                continue

        i += 1

    for v in results_dict.values():
        v["ex_conv"] = v["ex_conv"] / v["reps"]
        v["ex_div"] = v["ex_div"] / v["reps"]
        try:
            v["r_c"] = v["obs_conv"] / v["ex_conv"]
        except ZeroDivisionError:
            v["r_c"] = 0.
        try:
            v["r_d"] = v["obs_div"] / v["ex_div"]
        except ZeroDivisionError:
            v["r_d"] = 0.
        if v["r_c"] > 1:
            v["p_c"] = v["sims_c_gt"] / v["reps"]
        if v["r_c"] < 1:
            v["p_c"] = v["sims_c_lt"] / v["reps"]
        if v["r_d"] > 1:
            v["p_d"] = v["sims_d_gt"] / v["reps"]
        if v["r_d"] < 1:
            v["p_d"] = v["sims_d_lt"] / v["reps"]


    return results_dict


def run_sim_sf(tree: Node, alignment: dict, branch_combs: list[tuple[tuple]], model: str="JTT",
               rates: dict=None, use_anc: int=1, reps: int=100, write: bool=False) -> dict:
    """
    function to run and parse expect convergence and divergence by simulation, using site-specific
    frequencies. Frequencies will be calculated from the original alignment.
    """
    results_dict = {comb: {"ex_conv": 0, "ex_div": 0, "obs_conv": 0, "obs_div": 0, "sims_c_lt": 0,
                           "sims_c_gt": 0, "sims_d_lt": 0, "sims_d_gt": 0, "p_c": 0., "p_d": 0.,
                           "r_c": 0., "r_d": 0., "reps": 0} for comb in branch_combs}
    ex_br_dict = {br: [] for x in branch_combs for br in x}
    obs_br_dict = {br: [] for x in branch_combs for br in x}
    add_subs(obs_br_dict, alignment)
    for comb in branch_combs:
        try:
            obs_res = count_conv_subs(comb, obs_br_dict)
            results_dict[comb]["obs_conv"] = obs_res["CONV"]
            results_dict[comb]["obs_div"] = obs_res["DIV"]
        except ValueError:
            sys.stderr.write(f"skipping combination {comb} as no subs on at least one branch\n")
            continue
    i = 0
    while i < reps:
        sim_aln = sim_alignment_site_freqs(tree, alignment, model, use_anc, rates)
        sim_aln = sq.insert_gaps_by_seq(alignment, sim_aln)
        if write:
            sq.write_fasta(sim_aln, f"rep_{i}.pep.fa")
        # print(sq.get_fasta_str(sim_aln))
        add_subs(ex_br_dict, sim_aln)

        for comb in branch_combs:
            try:
                ex_res = count_conv_subs(comb, ex_br_dict)
                results_dict[comb]["ex_conv"] += ex_res["CONV"]
                results_dict[comb]["ex_div"] += ex_res["DIV"]
                if ex_res["CONV"] >= results_dict[comb]["obs_conv"]:
                    results_dict[comb]["sims_c_gt"] += 1
                if ex_res["CONV"] <= results_dict[comb]["obs_conv"]:
                    results_dict[comb]["sims_c_lt"] += 1
                if ex_res["DIV"] >= results_dict[comb]["obs_div"]:
                    results_dict[comb]["sims_d_gt"] += 1
                if ex_res["DIV"] <= results_dict[comb]["obs_div"]:
                    results_dict[comb]["sims_d_lt"] += 1
                results_dict[comb]["reps"] += 1
            except ValueError:
                sys.stderr.write(f"skipping combination {comb} as no subs on at least on branch\n")
                continue

        i += 1

    for v in results_dict.values():
        v["ex_conv"] = v["ex_conv"] / v["reps"]
        v["ex_div"] = v["ex_div"] / v["reps"]
        try:
            v["r_c"] = v["obs_conv"] / v["ex_conv"]
        except ZeroDivisionError:
            v["r_c"] = 0.
        try:
            v["r_d"] = v["obs_div"] / v["ex_div"]
        except ZeroDivisionError:
            v["r_d"] = 0.
        if v["r_c"] > 1:
            v["p_c"] = v["sims_c_gt"] / v["reps"]
        if v["r_c"] < 1:
            v["p_c"] = v["sims_c_lt"] / v["reps"]
        if v["r_d"] > 1:
            v["p_d"] = v["sims_d_gt"] / v["reps"]
        if v["r_d"] < 1:
            v["p_d"] = v["sims_d_lt"] / v["reps"]


    return results_dict


def run_sim_sf_indelible(tree: Node, alignment: dict, branch_combs: list[tuple[tuple]],
                         model: str="JTT", rates: dict=None, site_freqs: dict=None,
                         reps: int=100) -> dict:
    """
    function to run and parse expect convergence and divergence by simulation, using site-specific
    frequencies. Frequencies will be calculated from the original alignment if not supplied in
    site_freqs
    """
    results_dict = {comb: {"ex_conv": 0, "ex_div": 0, "obs_conv": 0, "obs_div": 0, "sims_c_lt": 0,
                           "sims_c_gt": 0, "sims_d_lt": 0, "sims_d_gt": 0, "p_c": 0., "p_d": 0.,
                           "r_c": 0., "r_d": 0., "reps": 0} for comb in branch_combs}
    ex_br_dict = {br: [] for x in branch_combs for br in x}
    obs_br_dict = {br: [] for x in branch_combs for br in x}
    add_subs(obs_br_dict, alignment)
    for comb in branch_combs:
        try:
            obs_res = count_conv_subs(comb, obs_br_dict)
            results_dict[comb]["obs_conv"] = obs_res["CONV"]
            results_dict[comb]["obs_div"] = obs_res["DIV"]
        except ValueError:
            sys.stderr.write(f"skipping combination {comb} as no subs on at least one branch\n")
            continue
    sim_reps = sim_sites_indelible(tree, alignment, model, rates, site_freqs, reps)
    for rep in sim_reps.values():
        add_subs(ex_br_dict, rep)

        for comb in branch_combs:
            try:
                ex_res = count_conv_subs(comb, ex_br_dict)
                results_dict[comb]["ex_conv"] += ex_res["CONV"]
                results_dict[comb]["ex_div"] += ex_res["DIV"]
                if ex_res["CONV"] >= results_dict[comb]["obs_conv"]:
                    results_dict[comb]["sims_c_gt"] += 1
                if ex_res["CONV"] <= results_dict[comb]["obs_conv"]:
                    results_dict[comb]["sims_c_lt"] += 1
                if ex_res["DIV"] >= results_dict[comb]["obs_div"]:
                    results_dict[comb]["sims_d_gt"] += 1
                if ex_res["DIV"] <= results_dict[comb]["obs_div"]:
                    results_dict[comb]["sims_d_lt"] += 1
                results_dict[comb]["reps"] += 1
            except ValueError:
                sys.stderr.write(f"skipping combination {comb} as no subs on at least on branch\n")
                continue

    for v in results_dict.values():
        v["ex_conv"] = v["ex_conv"] / v["reps"]
        v["ex_div"] = v["ex_div"] / v["reps"]
        try:
            v["r_c"] = v["obs_conv"] / v["ex_conv"]
        except ZeroDivisionError:
            v["r_c"] = 0.
        try:
            v["r_d"] = v["obs_div"] / v["ex_div"]
        except ZeroDivisionError:
            v["r_d"] = 0.
        if v["r_c"] > 1:
            v["p_c"] = v["sims_c_gt"] / v["reps"]
        if v["r_c"] < 1:
            v["p_c"] = v["sims_c_lt"] / v["reps"]
        if v["r_d"] > 1:
            v["p_d"] = v["sims_d_gt"] / v["reps"]
        if v["r_d"] < 1:
            v["p_d"] = v["sims_d_lt"] / v["reps"]


    return results_dict


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="tree to conduct simulations on, including node labels")
    parser.add_argument("alignment", help="alignment including ancestral sequences matching the \
                        internal node labels of tree to calculate observed counts")
    parser.add_argument("reps", help="number of simulation reps to conduct", type=int, default=100)
    parser.add_argument("branches", help="space-separated parent-daughter comparisons in format \
                        parent,daughter (at least two)", type=str, nargs="+")
    group = parser.add_mutually_exclusive_group(required = False)
    group.add_argument("-gf", "--gene_frequencies", help="use gene-wide empirical frequencies",
                        action="store_true")
    group.add_argument("-sf", "--site_frequencies", help="use site-specific empirical frequencies",
                        action="store_true")
    parser.add_argument("-sff", "--site_frequencies_file", help="if -sf, read frequencies from \
                        file (tsv, first column position, second column comma-separated \
                        frequencies)")
    parser.add_argument("-r", "--rates", help="per-site posterior mean rates from paml")
    parser.add_argument("-g", "--gamma", help="for simulations with gene-wide frequencies, also \
                        incorporate gamma-distributed rate heterogeneity. Two arguments, \
                        space-separated, ncat alpha", type=str, nargs=2)
    parser.add_argument("-i", "--indelible", help="use INDELible for the simulation under \
                        site-specific frequencies", action="store_true")
    parser.add_argument("-w", "--write", help="write individual replicate FASTAs. Note that \
                        INDELible", action="store_true")
    args = parser.parse_args()
    if args.site_frequencies_file:
        args.site_frequencies = True

    curroot = nwk.parse_from_file(args.tree)
    aln = dict(sq.parse_fasta(args.alignment))
    nodes = get_node_dict_from_tree(curroot)

    if args.rates:
        rates_dict = sq.parse_paml_rates(args.rates)
    else:
        rates_dict = None

    branches = [(x.split(",")[0], x.split(",")[1]) for x in args.branches]
    if len(branches) == 1:
        sys.stderr.write("More than one branch required\n")
        sys.exit(1)

    for branch in branches:
        if nodes[branch[0]] != nodes[branch[1]].parent:
            sys.stderr.write(f"branch {branch} is multi-branch lineage. Check.\n")

    if len(branches) == 2:
        br_combs = [tuple(branches)]
    else:
        br_combs = list(combinations(branches, 2))

    good_combs = get_good_branch_combs(br_combs, curroot)

    if args.site_frequencies:
        if args.site_frequencies_file:
            freqs = sq.parse_site_frequencies_file(args.site_frequencies_file)
        else:
            freqs = None
        if args.indelible:
            res = run_sim_sf_indelible(curroot, aln, good_combs, model="JTT", rates=rates_dict,
                                    site_freqs=freqs, reps=args.reps)
        else:
            res = run_sim_sf(curroot, aln, good_combs, model="JTT", rates=rates_dict, use_anc=1,
                             reps=args.reps, write=args.write)

    else:
        res = run_sim_gf(curroot, aln, good_combs, model="JTT", gamma=args.gamma, use_anc=1,
                         reps=args.reps, write=args.write)

    print("branches\tex_conv\tex_div\tobs_conv\tobs_div\tsims_c_lt\tsims_c_gt\tsims_d_lt\t"
          "sims_d_gt\tp_c\tp_d\tr_c\tr_d\treps")
    for br_comb, nums in res.items():
        branches = " ".join([",".join(b) for b in br_comb])
        print(f"{branches}\t" + '\t'.join([f"{x}" for x in nums.values()]))
