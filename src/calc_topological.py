#! /usr/bin/python3

import sys, math, subprocess
import pandas as pd

# requires iqtree, but you should be able to substitute raxml just as well

def tp_to_dict(inf):
    """converts TREE-PUZZLE formatted site loglikelihoods into a dictionary for use with pandas. Returns a pandas df."""
    with open(inf,"r") as tp:
        header = tp.readline()
        nsites = int(header.strip().split()[1])
        sites = range(1,nsites+1,1)
        ll = [float(x) for x in tp.readline().strip().split()[1:]]
        tmp_dict = {"sites":sites,"logl":ll}
        df = pd.DataFrame(tmp_dict)
        return df

def optimize_model(aln,intree,pre,model=None):
    """optimizes branch lengths and model params for the specified alignment and topology for the best-fitting aa substitution model."""
    if model is not None:
        cmd = ["iqtree","-s",aln,"-m",model,"-te",intree,"-nt","1","-pre",pre,"-wsl","-redo"]
    else:
        cmd = ["iqtree","-s",aln,"-m","MFP","-msub","nuclear","-mrate","E,G","-te",intree,"-nt","1","-pre",pre,"-wsl","-redo"] 
    p = subprocess.Popen(cmd,shell=False)
    p.wait()
    with open(pre+".iqtree","r") as mfile:
        for line in mfile:
            line = line.strip()
            if "Model of substitution" in line:
                mstring = line.split(":")[1].strip()
    if model is None:
        return mstring

def prob_ap(x,y): # from https://gitlab.in2p3.fr/pveber/reviewphiltrans/blob/master/lib/scripts/calc_topological.py
    if x and y:
        return math.exp(x-y)/(math.exp(x-y)+1)
    else:
        return 0

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("usage: python "+sys.argv[0]+" aln tree conv_tree outf")
        sys.exit()

    m = optimize_model(sys.argv[1],sys.argv[2],"non_conv",None) # non-convergent topology
    optimize_model(sys.argv[1],sys.argv[3],"conv",m) # convergent topology

    df_non_conv = tp_to_dict("non_conv.sitelh") 
    """in general bits below are based on https://gitlab.in2p3.fr/pveber/reviewphiltrans/blob/master/lib/scripts/calc_topological.py"""
    df_conv = tp_to_dict("conv.sitelh")

    df_final = pd.merge(df_non_conv[["sites","logl"]], df_conv[["sites","logl"]],
                 on = "sites",
                 suffixes=['_non_conv', '_conv'])

    df_final["topological"] = list(map(prob_ap,  df_final["logl_conv"],df_final["logl_non_conv"]))
    df_final = df_final[["sites","topological"]]
    df_final.to_csv(sys.argv[4],sep=",",index=False)



