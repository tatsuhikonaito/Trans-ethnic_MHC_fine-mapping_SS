#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import subprocess
import numpy as np
import pandas as pd
from scipy import stats


# R2 threshold for variant filtering
thr = 0.7


def ss_conditional_analysis(args):

    dishfile_list = pd.read_csv(args.dishfile_list, sep="\t", header=None, names=["dishfilename", "n_case", "n_control"])
    dishfile_list["n_mean"] = stats.hmean((np.array(dishfile_list["n_case"]), np.array(dishfile_list["n_control"])))
    ref_list = pd.read_csv(args.ref_list, sep="\t", header=None, names=["ref_filename"])
    n_study = len(dishfile_list)

    # Integrate summary files
    for i in range(n_study):
        fn = dishfile_list.loc[i, "dishfilename"]
        dishfile = pd.read_csv(fn, sep="\t", header=0)
        dishfile = dishfile[dishfile.r2pred>=thr]

        ref = pd.read_csv(ref_list.loc[i, "ref_filename"]+".bim", sep="\t", header=None, names=['chr', 'id', 'dist', 'pos', 'a1', 'a2'])
        subprocess.call("plink --bfile {0} --freq --out {0}".format(ref_list.loc[i, "ref_filename"]), shell=True)
        freq = pd.read_csv(ref_list.loc[i, "ref_filename"]+".frq", sep="\t|\s+", header=0, index_col=1, engine="python")

        dishfile["N"] = dishfile_list.loc[i, "n_mean"]
        dishfile["maf"] = freq.loc[dishfile.Marker_id, "MAF"]
        sd = np.std([1]*dishfile_list.loc[i, "n_case"] + [0]*dishfile_list.loc[i, "n_control"])
        c = 1/sd
        dishfile["BETA"] = c*(dishfile["Imputed_Z"]/((2*dishfile["maf"]*(1-dishfile["maf"])*(dishfile["N"]+dishfile["Imputed_Z"]**2))**(1/2)))
        dishfile["SE"] = dishfile["BETA"]/dishfile["Imputed_Z"]

        dishfile.to_csv("{0}.ma".format(fn), sep="\t", columns=["Marker_id", "Effect_allele", "Non_effect_allele", "maf", "BETA", "SE", "imputed_P", "N"], \
                       header=["SNP", "A1", "A2", "freq", "b", "se", "p", "N"], index=False)
        subprocess.call("gcta64 --bfile {0} --cojo-file {1}.ma --cojo-cond {2} --out {1}".format(ref_list.loc[i, "ref_filename"], fn, args.allele_list), shell=True)
        

def main():
    parser = argparse.ArgumentParser(description='Perform sample size-based meta-analysis of Z-score.')
    parser.add_argument('--dishfile-list', required=True, help='List of dish output filenames', dest='dishfile_list')
    parser.add_argument('--ref-list', required=True, help='List of HLA imputation reference filenames', dest='ref_list')
    parser.add_argument('--allele-list', required=True, help='List of alleles to be condition in', dest='allele_list')
    args = parser.parse_args()

    ss_conditional_analysis(args)


if __name__ == '__main__':
    main()
