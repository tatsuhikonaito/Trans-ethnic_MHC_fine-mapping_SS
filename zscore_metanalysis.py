#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import numpy as np
import pandas as pd


# R2 threshold for variant filtering
thr = 0.7


def zscore_metaanalysis(args):

    file_list = pd.read_csv(args.file_list, sep="\t", header=None, names=["finename", "n"])
    n_study = len(file_list)
    n_total = np.sum(file_list["n"])

    # Integrate summary files
    for i in range(n_study):
        fn = file_list.loc[i, "filename"]
        imputed = pd.read_csv(fn, sep="\t", header=0)
        imputed = imputed[imputed.r2pred>=thr]
        if i == 0:
            all_allele = 
            result = pd.DataFrame(index=imputed.index, columns=["z_{}".format(i+1) for i in range(n_study)]+["z_meta", "p_meta"])
            result["z_1"] = imputed["Imputed_Z"]
        else:
            result = result.loc[result.index[result.index.isin(imputed.index)]]
            result["z_{}".format(i+1)] = imputed.loc[result.index, "Imputed_Z"]

    # Z-score meta-analysis
    for i, r in result.iterrows():
        result.loc[i, "z_meta"] = np.dot(r["z_{}".format(i+1) for i in range(n_study)], (file_list["n"]/n_total)**(1/2))
        result.loc[i, "p_meta"] = 2*norm.sf(np.abs(result.loc[i, "z_meta"]))

    result.to_csv("metaanalysis_result.txt", index=True, header=True, sep="\t")
        

def main():
    parser = argparse.ArgumentParser(description='Perform sample size-based meta-analysis of Z-score.')
    parser.add_argument('--file-list', required=True, help='List of summary statistics filenames', dest='file_list')

    args = parser.parse_args()
    zscore_metaanalysis(args)


if __name__ == '__main__':
    main()
