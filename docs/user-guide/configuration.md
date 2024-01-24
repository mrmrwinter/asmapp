# Configuration

In order to run the workflow, the user needs to set up a configuration file, `config.yaml`. Open this in a text editor and enter information according to the annotations. 

Next, place the assembly you want to be appraised, in `.fasta` format, in the `data/assemblies/` directory, and place the reads used to generate it, in `.fastq.gz` format, in the `data/reads/` directory. This will be streamlined in a future release.

Singularity is required for several tages of the workflow. This can be installed from [here](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

Contaminant detection with __Blobtools__ requires installation of the entire NCBI nt database. Instructions for this can be found [here](https://ftp.ncbi.nlm.nih.gov/blast/db/), with instructions found [here](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html). 

In order to compare to a reference, and detect mitochondrial contigs, ASMAPP needs sequences to compare against.
Download and place a reference nuclear genome sequence, and a reference mitochondrial sequence, in the `data/assemblies/` directory. If these analyses are not required, ignore this step, and add `-k` to the main snakemake command when running the workflow.
