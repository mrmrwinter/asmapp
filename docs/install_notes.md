# Notes on installation

## MINICONDA

This workflow requires miniconda. Miniconda can be downloaded here: https://docs.conda.io/en/latest/miniconda.html

<br>

## BLOBPLOTS

Blobplots requires some databases installed locally in order to perform a taxonomy lookup.

I will add a setting in the config to path to a locally installed blast database, as well as instructions on how to install one, here, but I'll also add a rule to install the blast database if none is present.

<br>

### BLAST database install instructions

setup requires that the ncbi database is installed and configured in the config files

NCBI nt database can be installed using these instructions: https://www.ncbi.nlm.nih.gov/books/NBK569850/

Also needs taxdump downloading and unpacking in the same folder as the nt database
To do this, navigate to the directory housing the nt database and run the following: `wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz`
And then extract it into the same directory

<br>

### QUAST installation

QUAST needs to be downloaded for the workflow to use it.
This can be done with the following line of code (make sure to run it somewhere convenient):
`git clone https://github.com/ablab/quast.git`

<br>

### Configuring the workflow

Now go to config_local.yaml and edit it acording to your needs and paths

<br>

### Running the workflow

To run the workflow navigate at command line to the asmapp directory, and enter the following code:
`snakemake -cX --configfile config_local.yaml --use-singularity -k`

