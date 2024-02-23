#! /usr/bin/python3


import os
import re
import sys
import argparse
import subprocess
import numpy as np


MODELS = ["Blosum62", "cpREV", "Dayhoff", "DCMut", "DEN", "FLU",
          "HIVb", "HIVw", "JTT", "JTT-DCMut", "LG", "mtART",
          "mtMAM", "mtREV", "mtZOA", "PMB", "rtREV", "stmtREV",
          "VT", "WAG", "LG4M", "LG4X", "PROTGTR"]


def optimise(tree, aln, model, outgroup=None):
    if model.split("+")[0] not in MODELS:
        # check if valid model
        sys.stderr.write("please specify a valid AA model!")
        sys.exit()
    if not set(model.split("+")[1:]).intersection(set(["F", "FO", "G"])):
        # check if valid suffix
        sys.stderr.write("only empirical frequencies (+F) \
                         and/or gamma-distributed rates (+G) \
                         allowed.\n")
        sys.exit()
    if outgroup is not None:
        # cmd with OG (drawing option only, produces rooted tree)
        cmd = ["raxml-ng",
               "--evaluate",
               "--tree",
               tree,
               "--msa",
               aln,
               "--model",
               model,
               "--prefix",
               model + "_tree",
               "--outgroup",
               ",".join(outgroup),
               "--lh-epsilon",
               "0.001"
               ]
    else:
        # no OG
        cmd = ["raxml-ng",
               "--evaluate",
               "--tree",
               tree,
               "--msa",
               aln,
               "--model",
               model,
               "--prefix",
               model + "_tree",
               "--lh-epsilon",
               "0.001"
               ]
    p = subprocess.Popen(cmd, shell=False,
                         stdout=subprocess.DEVNULL,
                         stderr=subprocess.DEVNULL)
    p.wait()


def get_optimised_freqs(tree: str, aln: str, model: str) -> np.array:
    if model.split("+")[0] not in MODELS:
        # check if valid model
        sys.stderr.write("please specify a valid AA model!\n")
        sys.exit()
    if "FO" not in model.split("+"):
        # check that we are calling for optimised freqs
        sys.stderr.write("this function is for optimised frequencies\n")
        sys.exit()
    cmd = ["raxml-ng",
           "--evaluate",
           "--tree",
           tree,
           "--msa",
           aln,
           "--model",
           model,
           "--opt-branches",
           "off",
           "--force",
           "model_lh_impr"]
    print(subprocess.list2cmdline(cmd))
    p = subprocess.Popen(cmd, shell=False,
                         stdout=subprocess.DEVNULL,
                         stderr=subprocess.DEVNULL)
    p.wait()
    with open(f"{aln}.raxml.bestModel", "r") as modelf:
        line = modelf.readline()
        line = line.strip().split(",")[0]
        freqs = line.split("{")[1].rstrip("}")
        freqs = [float(x) for x in freqs.split("/")]
        freqs = np.array(freqs)
        freqs = freqs / sum(freqs)
    
    # tidy up
    pat = re.compile(fr"{aln}.raxml\..*")
    for f in os.listdir():
        if re.search(pat, f):
            os.remove(f)
        elif f == aln:
            os.remove(f)
        
    return freqs




def root_and_clean(tree, outgroup):
    rrCmd = ["pxrr",
             "-t",
             tree,
             "-g",
             ",".join(outgroup),
             "-o",
             tree + ".rr"]
    print(rrCmd)
    cltrCmd = ["pxcltr",
               "-t",
               tree + ".rr",
               "-o",
               tree + ".rr.cltr"]
    print(cltrCmd)
    # incase tree has annotations like from GeneRax
    subprocess.Popen(rrCmd, shell=False,
                     stdout=subprocess.DEVNULL,
                     stderr=subprocess.DEVNULL).wait()
    subprocess.Popen(cltrCmd, shell=False,
                     stdout=subprocess.DEVNULL,
                     stderr=subprocess.DEVNULL).wait()


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="nwk tree file to optimise model \
                        and brlen")
    parser.add_argument("aln", help="FASTA alignment")
    parser.add_argument("model", help="Model string, raxml-ng")
    parser.add_argument("-og", "--outgroup", help="comma-separated list \
                        of OG taxa to root on", type=str)
    args = parser.parse_args()
    if args.outgroup is not None:
        OG = [x for x in args.outgroup.split(",")]
    else:
        OG = None

    optimise(args.tree, args.aln, args.model, OG)
    if args.outgroup is not None:
        root_and_clean(args.model + "_tree.raxml.bestTree", OG)
