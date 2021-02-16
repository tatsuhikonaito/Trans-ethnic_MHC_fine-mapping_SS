#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import numpy as np
import pandas as pd
from scipy import stats


# R2 threshold for variant filtering
thr = 0.7


def zscore_metaanalysis(args):

    dishfile_list = pd.read_csv(args.dishfile_list, sep="\t", header=None, names=["dishfilename", "n_case", "n_control"])
    dishfile_list["n_mean"] = stats.hmean((np.array(dishfile_list["n_case"]), np.array(dishfile_list["n_control"])))
    n_study = len(dishfile_list)
    n_total = np.sum(dishfile_list["n_mean"])

    # Integrate summary files
    for i in range(n_study):
        dishfile = pd.read_csv(dishfile_list.loc[i, "dishfilename"], sep="\t", header=0, index_col=0)
        dishfile = dishfile[dishfile.r2pred>=thr]
        if i == 0:
            result = pd.DataFrame(index=dishfile.index, columns=["z_{}".format(i+1) for i in range(n_study)]+["z_meta", "p_meta"])
            result["z_1"] = dishfile["Imputed_Z"]
        else:
            result = result.loc[result.index[result.index.isin(dishfile.index)]]
            result["z_{}".format(i+1)] = dishfile.loc[result.index, "Imputed_Z"]

    # Z-score meta-analysis
    for i in result.index:
        result.loc[i, "z_meta"] = np.dot(result.loc[i, ["z_{}".format(i+1) for i in range(n_study)]], (dishfile_list["n_mean"]/n_total)**(1/2))
        result.loc[i, "p_meta"] = 2*stats.norm.sf(np.abs(result.loc[i, "z_meta"]))

    result.to_csv("metaanalysis_result.txt", index=True, header=True, sep="\t")
        

def main():
    parser = argparse.ArgumentParser(description='Perform sample size-based meta-analysis of Z-score.')
    parser.add_argument('--dishfile-list', required=True, help='List of dish output filenames', dest='dishfile_list')
    args = parser.parse_args()

    zscore_metaanalysis(args)


if __name__ == '__main__':
    main()
