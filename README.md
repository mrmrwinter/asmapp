# ASMAPP

## This workflow is for appraising genome assemblies. 


### This workflow includes full assembly appraisal with the following methods:  
-Read mapping   
-Assembly appraisal statistics  
-Mitochondrial contig flagging/removal  
-Contaminant identification   
-Homeologous contig detection   
-Statistics on homeolog differences  
-Variant calling   
-Mapping stat generation   
-Ploidy estimations   
-Coverage analytics   

More detailed instructions about how to run the workflow will be included in a docs file at a later date.

<br>

## Installation and running instructions, for now:
  
1, Install conda
  
2, Clone this repository and navigate to its top level directory
  
3, Run `conda env create -f envs/asmapp.yaml`
  
4, Activate the environment - `conda activate asmapp`
  
5, Configure the `config.yaml`
  
6, Run the workflow with the follownig command; `snakemake -cX -use-singularity`, whre X is the same number of cores desired.

Full dependency installation instructions can be found at `docs/user-guide/installation.md`
  

