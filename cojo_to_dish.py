#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd


def cojo_to_dish(args):

    cojofile_list = pd.read_csv(args.cojofile_list, sep="\t", header=None, names=["cojofilename"])
    ref_list = pd.read_csv(args.ref_list, sep="\t", header=None, names=["ref_filename"])

    for i in cojofile_list.index:
        cojofile = pd.read_csv(cojofile_list.loc[i, "cojofilename"], sep="\t", header=0)
        ref = pd.read_csv(ref_list.loc[i, "ref_filename"]+".bim", sep="\t", header=None, names=['chr', 'id', 'dist', 'pos', 'a1', 'a2'])

        cojofile["z"] = cojofile["bC"]/cojofile["bC_se"]
        cojofile["r2pred"] = 1  # R2preds of DISH of all variants are set to 1 for convenience.
        cojofile["a1"] = ref.loc[cojofile.index, "a1"]
        cojofile["a2"] = ref.loc[cojofile.index, "a2"]

        dishfile = cojofile[["SNP", "bp", "a1", "a2", "z", "r2pred", "pC"]]
        dishfile.columns = ["Marker_id", "Marker_pos", "Effect_allele", "Non_effect_allele", "Imputed_Z", "r2pred", "imputed_P"]

        dishfile.to_csv("{0}.dish.txt".format(cojofile_list.loc[i, "cojofilename"].rstrip(".cma.cojo")), sep="\t", header=True, index=False)
        

def main():
    parser = argparse.ArgumentParser(description='Perform sample size-based meta-analysis of Z-score.')
    parser.add_argument('--cojofile-list', required=True, help='List of cojo output filenames', dest='cojofile_list')
    parser.add_argument('--ref-list', required=True, help='List of HLA imputation reference filenames', dest='ref_list')
    args = parser.parse_args()

    cojo_to_dish(args)


if __name__ == '__main__':
    main()
