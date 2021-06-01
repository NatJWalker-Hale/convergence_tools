#! /usr/bin/python3

import sys
import argparse
import subprocess
from parse_fasta import parse_fasta


def rename_seqs(seqDict1, seqDict2):
    names = {}
    rnSeqDict1 = {}
    rnSeqDict2 = {}
    count = 0
    for k, _ in seqDict1.items():
        rnSeqDict1["1_S"+str(count)] = seqDict1[k]
        names["1_S"+str(count)] = k
        count += 1
    count = 0
    for k, _ in seqDict2.items():
        rnSeqDict2["2_S"+str(count)] = seqDict2[k]
        names["2_S"+str(count)] = k
        count += 1
    return rnSeqDict1, rnSeqDict2, names


def get_site_dict(seqDict):
    siteDict = {}
    nSite = len([v for v in seqDict.values()][0])
    s = 0
    while s < nSite:
        states = [x[s] for x in seqDict.values()]
        siteDict[s] = states
        s += 1
    return siteDict


def align(seqDict1, seqDict2):
    with open("tmp1", "w") as outf:
        for k, v in seqDict1.items():
            outf.write(">"+str(k)+"\n")
            outf.write(v+"\n")
    with open("tmp2", "w") as outf:
        for k, v in seqDict2.items():
            outf.write(">"+str(k)+"\n")
            outf.write(v+"\n")
    cmd = ["muscle", "-profile", "-in1", "tmp1", "-in2", "tmp2", "-out",
           "tmp.aln"]
    subprocess.run(cmd)


def check_aligned(seqDict):
    lens = []
    for v in seqDict.values():
        lens.append(len(v))
    if len(set(lens)) > 1:
        return False
    else:
        return True


def split_out_alignment(outAlnDict):
    outAln1 = {}
    outAln2 = {}
    for k, v in outAlnDict.items():
        if k.startswith("1_"):
            outAln1[k] = v
        elif k.startswith("2_"):
            outAln2[k] = v
    return outAln1, outAln2


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("aln1", help="first alignment")
    parser.add_argument("aln2", help="second alignment")
    args = parser.parse_args()

    seq1 = dict([x for x in parse_fasta(args.aln1)])
    if not check_aligned(seq1):
        sys.stderr.write("sequences in aln1 not aligned\n")
        sys.exit()
    seq2 = dict([x for x in parse_fasta(args.aln2)])
    if not check_aligned(seq2):
        sys.stderr.write("sequences in aln2 not aligned\n")
        sys.exit()

    rnSeq1, rnSeq2, nameCorres = rename_seqs(seq1, seq2)
    align(rnSeq1, rnSeq2)
    out = dict([x for x in parse_fasta("tmp.aln")])

    out1, out2 = split_out_alignment(out)
    len1 = len([
                x for x in rnSeq1.values()
                ][0])
    len2 = len([
                x for x in rnSeq2.values()
                ][0])

    """we can exploit the fact that MUSCLE profile
    will only either keep the length of the shorter
    alignment, or insert whole columns of gaps"""

    startDiff = 0
    corres = {}
    if len1 <= len2:
        outCols1 = get_site_dict(out1)
        s = 0
        first = False
        for k, v in outCols1.items():
            chars = set(v)
            if len(chars) == 1 and "-" in chars:
                # whole column of gaps
                if first:
                    startDiff += 1
                    corres[s] = k - startDiff
                else:
                    startDiff += 1
            else:
                first = True
                corres[s] = k - startDiff
            s += 1
    for k, v in corres.items():
        print(str(k+1) + "\t" + str(v+1))
