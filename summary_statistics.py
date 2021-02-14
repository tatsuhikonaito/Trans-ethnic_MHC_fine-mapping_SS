#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import subprocess
import numpy as np
import pandas as pd


# R2 threshold for variant filtering
thr = 0.7


def ss_conditional_analysis(args):

    file_list = pd.read_csv(args.file_list, sep="\t", header=None, names=["finename", "n"])
    ref_list = pd.read_csv(args.ref_list, sep="\t", header=None, names=["ref_filename"])
    n_study = len(file_list)
    n_total = np.sum(file_list["n"])

    # Integrate summary files
    for i in range(n_study):
        fn = file_list.loc[i, "filename"]
        imputed = pd.read_csv(fn, sep="\t", header=0)
        imputed = imputed[imputed.r2pred>=thr]

        ref = pd.read_csv(ref_list.loc[i, "ref_filename"]+".bim", sep="\t", header=None, names=['chr', 'id', 'dist', 'pos', 'a1', 'a2'])
        subprocess.call("plink --bfile {0} --freq --out {}".format(ref_list.loc[i, "ref_filename"]), shell=True)
        freq = pd.read_csv(ref_list.loc[i, "ref_filename"]+".frq", sep="\t|\s+", header=0, engine="python")

        imputed["N"] = file_list.loc[i, "n"]
        imputed["maf"] = freq.loc[imputed.Marker_id, "MAF"]
        imputed["BETA"] = c*(imputed["Imputed_Z"]/((2*imputed["maf"]*(1-imputed["maf"])*(n+imputed["Imputed_Z"]**2))**(1/2)))
        imputed["SE"] = imputed["BETA"]/imputed["imputed_Z"]

        imputed.to_csv("{0}.ma".format(fn), sep="\t", columns=["Marker_id", "Effect_allele", "Non_effect_allele", "maf", "BETA", "SE", "imputed_P", "N"], \
                       header=["SNP", "A1", "A2", "freq", "b", "se", "p", "N"], index=False)
        subprocess.call("gcta64 --bfile {0} --cojo-file {1}.ma --cojo-cond {2} --out {0}.{1}".format(ref_list.loc[i], fn, args.allele_list), shell=True)
        

def main():
    parser = argparse.ArgumentParser(description='Perform sample size-based meta-analysis of Z-score.')
    parser.add_argument('--file-list', required=True, help='List of summary statistics filenames', dest='file_list')
    parser.add_argument('--ref-list', required=True, help='List of HLA imputation reference filenames', dest='ref_list')
    parser.add_argument('--allle-list', required=True, help='List of alleles to be condition in', dest='allele_list')
    args = parser.parse_args()

    zscore_metaanalysis(args)


if __name__ == '__main__':
    main()
