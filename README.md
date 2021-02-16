# Trans-ethnic fine-mapping of the MHC region using GWAS summary statistics
This repository is a storage of source codes used in our study of trans-ethnic MHC fine-mapping on Parkinson's disease risk.

## Publication/Citation

The study is described in the following paper. 

- N.A.

Please cite this paper if you use any material in this repository.

## Requirements

- Python 3.x (3.7.4) with the following modules.
  - Numpy (1.17.2)
  
  - Pandas (0.25.1)
  
  - Scipy (1.3.1)
  
  - Argparse (1.4.0)
  
    Our scripts were tested on the versions in parentheses, so we do not guarantee that it will work on different versions.
  
- R

- Plink [1]

- DISH [2]

- GCTA-COJO [3]

## Installation

Just clone this repository as folllows.

```
git clone https://github.com/tatsuhikonaito/Trans-ethnic_MHC_finemapping_SS
cd ./Trans-ethnic_MHC_finemapping_SS
```

## Process

### 0. Preparation

The access to the GWAS summary statistics and HLA imputation reference panels (Plink binary format) used in our study is described in the paper.

### 1. Fine-mapping in the MHC region using summary statistics

Run `DISH` to perform MHC fine-mapping for post-QC summary statistics data of each study.

```
$ Rscript DISH.r <summary filename> T hg19 EUR 0.005 T <output filename> 0.05 
```

##### Inputs

- \<summary filename>

  Filename of GWAS summary statistics to impute

##### Outputs

- \<output filename>

  Filename of DISH output of HLA imputation.

### 2. Sample size-based meta-analysis of Z-score

Peform sample size-based meta-analysis of Z-score using imputed summary statistics files as follows. 

```
$ python zcore_metaanalysis.py --dishfile-list DISHFILE_LIST.txt
```

##### Inputs

- DISHFILE_LIST.txt

  First and second columns are the DISH-output filename and number of cases and controls of each study. 

  e.g.)

  study_1.dish.txt	1000	10000

  study_2.dish.txt	2000	10000

##### Outputs

- metaanalysis_result.txt

  Second to last and last columns are Z-score and P-value of meta-analysis, following Z-scores of individual studies.

### 3. Summary statistics conditional analysis

Peform sample size-based meta-analysis of Z-score using imputed summary statistics files as follows. 

```
$ python ss_conditional_analysis.py --dishfile-list DISHFILE_LIST.txt --ref-list REF_LIST.txt --allele-list ALLELE_LIST.txt
```

##### Inputs

- DISHFILE_LIST.txt

  Described above.

- REF_LIST.txt

  The filenames (prefix) of HLA imputation reference panels corresponding to individual studies in the same order as DISHFILE_LIST.txt.

  e.g.)

  T1DGC_REF

  PAN-Asian_REF

- ALLELE_LIST.txt

  First column is the name of alleles to be conditioned in.

  e.g.)

  AA_DRB1_13_32660109_R
  AA_DRB1_13_32660109_H

#### Outputs

- GCTA-COJO output files

  Prefix of filename of each study + ".cma.cojo".

### Conversion of file formats

Convert GCTA-COJO output file format to DISH output file format.

This enables iterative process of Z-score meta-analysis and conditional analysis.

```
$ python cojo_to_dish.py --cojofile-list COJOFILE_LIST.txt --ref-list REF_LIST.txt
```

##### Inputs

- COJOFILE_LIST.txt

  The filenames of HLA imputation reference panels corresponding to individual studies in the same order as DISHFILE_LIST.txt.

  e.g.)

  study_1.cma.cojo

  study_2.cma.cojo

- REF_LIST.txt

  Described above.

#### Outputs

- DISH output files

  Prefix of filename of each study + ".dish.txt".

## References

[1] Purcell S, Neale B, Todd-Brown K, et al. PLINK: a tool set for whole-genome association and population-based linkage analyses. Am. J. Hum. Genet. 2007;81(3):559–575.

[2] Lim J, Bae S-C, Kim K. Understanding HLA associations from SNP summary association statistics. Sci. Rep. 2019;9(1):1337.

[3] Yang J, Ferreira T, Morris AP, et al. Conditional and joint multiple-SNP analysis of GWAS summary statistics identifies additional variants influencing complex traits. Nat. Genet. 2012;44(4):369–375.

## Contact

For any question, you can contact Tatsuhiko Naito ([tnaito@sg.med.osaka-u.ac.jp](mailto:tnaito@sg.med.osaka-u.ac.jp))
