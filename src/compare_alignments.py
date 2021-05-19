#! /usr/bin/python3

import sys
import argparse
from parse_fasta import parse_fasta


def rename_seqs(seqDict1, seqDict2):
    count = 0
    for k, _ in seqDict1.items():
        seqDict1["1_S"+str(count)] = seqDict1.pop(k)
    for k, _ in seqDict2.items():
        seqDict2["2_S"+str(count)] = seqDict2.pop(k)


