#!/usr/bin/env python3


"""
bootstrap expected number of convergences and divergences for a given alignment and tree using 
simulation
"""

import sys
import argparse
import subprocess
from dataclasses import dataclass, field
import numpy as np
import sequence as sq
import newick as nwk
from phylo import Node
from sim_site_specific_freqs import sim_alignment_gene_freqs, sim_alignment_site_freqs
from calc_expected_conv import get_good_branch_combs
from count_conv_subs_from_asr import read_subs, count_conv_subs



@dataclass
class BranchComb:
    """
    simple class to hold simulation results for single branch combo
    """
    branches: tuple[tuple] = field(default_factory=tuple)
    conv: int=0
    div: int=0


@dataclass
class Simulation:
    """
    simple class to hold results of single simulation
    """
    branch_combs: list[BranchComb] = field(default_factory=list)
    alnlength: int=0
    sites: int=0
    

@dataclass
class SimulationResult:
    """
    simple class to hold results of n iterations of simulation
    """
    reps: int=0
    branch_combs: list[BranchComb] = field(default_factory=list)


def run_sim_gf(tree: Node, alignment: dict, branch_combs: tuple[tuple], model: str="JTT",
               gamma: tuple[int, float]=None, frequencies: np.array=None, use_anc: int=1,
               reps: int=100) -> SimulationResult:
    """
    function to run and parse expect convergence and divergence by simulation, using gene-wide
    frequencies. Frequencies can be specified or if not will be calculated from the input alignment
    """
    results = SimulationResult()
    i = 0
    while i < reps:
        sim_aln = sim_alignment_gene_freqs(tree, alignment, model, gamma, use_anc)
        sim_aln = sq.insert_gaps_by_seq(alignment, sim_aln)
        
        for comb in branch_combs:
            brc_result = BranchComb()
            brc_result.branches = comb
