# Quick start
  
1, Clone this repository and navigate to its top level directory
> `git clone https://github.com/mrmrwinter/asmapp.git`  
> `cd asmapp`
  
2, Create and activate the environment with conda
> `conda env create -f envs/asmapp.yaml`  
> `conda activate asmapp`

3, Download databases and packages. See `docs/user-guide/installation.md`

4, Insert inputs
- Place the assembly fasta and the reads into the respective folders, data/assemblies and data/reads

5, Configure the `config.yaml`
  
6, Run the workflow, where X is the number of cores to use, then generate a report  
> `snakemake -cX --use-singularity --use-conda`  
> `snakemake --report report.html`


  
