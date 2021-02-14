# Trans-ethnic fine-mapping of the major histocompatibility complex region using GWAS summary statistics
This repository is a storage of source codes used in our study of trans-ethnic fine-mapping of MHC region on Parkinson's disease risk.

## Publication/Citation

The study is described in the following paper. 

- N.A.

Please cite this paper if you use any material in this repository.

## Requirements

- Python 3.x (3.7.4) with the following modules (our scripts were tested on the versions in parentheses, so we do not guarantee that it will work on different versions.)
  - Numpy (1.17.2)
  - Pandas (0.25.1)
  - Scipy (1.3.1)
  - Argparse (1.4.0)
- R
- DISH [1]
- GCTA-COJO [2]

## Installation

Just clone this repository as folllows.

```
git clone https://github.com/tatsuhikonaito/DEEP-HLA
cd ./DEEP-HLA
```

## Process

### 1. Fine-mapping in the MHC region using summary statistics

Run `DISH` to perform MHC fine-mapping for post-QC summary statistics data of each study.

```
$ Rscript DISH.r <summary filename> T hg19 EUR 0.005 T <output filename> 0.05 
```


### 2. Sample size-based meta-analysis of Z-score

Peform sample size-based meta-analysis of Z-score using imputed summary statistics files as follows. 

```
$ python zcore_metaanalysis.py --file-list FILE_LIST.txt
```

##### Inputs

- FILE_LIST.txt

  First and second columns are the DISH-result filename and number of samples (harmonic mean) of each study. 

  e.g.)

  study_1.txt	1000

  study_2.txt	2000

##### Outputs

- metaanalysis_result.txt

  Second to last and last columns are Z-score and P-value of meta-analysis, following Z-scores of individual studies.

### 3. Summary statistics conditional analysis

Peform sample size-based meta-analysis of Z-score using imputed summary statistics files as follows. 

```
$ python ss_conditional_analysis.py --file-list FILE_LIST.txt --ref-list REF_LIST.txt --allele-list ALLELE_LIST.txt
```

- ##### Inputs

  - FILE_LIST.txt

    Described above.

  - REF_LIST.txt

    The filenames of HLA imputation reference panels corresponding to individual studies in the same order as FILE_LIST.txt.

    e.g.)

    T1DGC_REF

    PAN-Asian_REF

  - ALLELE_LIST.txt

    First column is the name of alleles HLA imputation reference panel corresponding to each study.

    e.g.)

    AA_DRB1_13_32660109_R
    AA_DRB1_13_32660109_H

- **Outputs:** GCTA-COJO result files.

### Others. Conversion of file formats

Peform sample size-b

## Reference

[1] Lim J, Bae S-C, Kim K. Understanding HLA associations from SNP summary association statistics. Sci. Rep. 2019;9(1):1337.

[2] Yang J, Ferreira T, Morris AP, et al. Conditional and joint multiple-SNP analysis of GWAS summary statistics identifies additional variants influencing complex traits. Nat. Genet. 2012;44(4):369–375.

## Contact

For any question, you can contact Tatsuhiko Naito ([tnaito@sg.med.osaka-u.ac.jp](mailto:tnaito@sg.med.osaka-u.ac.jp))
