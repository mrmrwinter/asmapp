
# Installation

## Downloading the workflow
To download the workflow, navigate to the directory you want the ASMAPP folder to be and run the following at the command line:
`git clone https://github.com/mrmrwinter/asmapp.git`
If you are unfamiliar with Git, or do not have it installed, download it [here](https://git-scm.com/downloads). 

## Dependencies
This workflow requires miniconda. Miniconda can be downloaded [here](https://docs.conda.io/en/latest/miniconda.html).
The majority of the dependencies required by this workflow are installed and managed through Miniconda. Packages and databases that require additional installation will be described below.

To install and activate the ASMAPP conda environment, navigate to the ASMAPP directory and run the following command:  
`conda create -f envs/asmapp.yaml -n asmapp`  

Following environment creation, run:  
`conda activate asmapp`


## Databases
 Some packages in the workflow require that the NCBI nt database is installed locally. This can be installed from the NCBI website [here](https://www.ncbi.nlm.nih.gov/books/NBK569850/). It is also necessary to donwload the NCBI `taxdb.tar.gz` from this [ftp](ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz) and extract it in the same folder as the nt database.

