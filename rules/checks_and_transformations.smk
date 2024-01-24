# INPUT CHECKS AND INITIAL DATA TRANSFORMATIONS

# Check for presence of input assembly
rule input_assembly:
    output:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",

# Check for presenc of input reads
rule input_reads:
    output:
        zipped_fastq_reads = "data/reads/" + config["reads"] + ".fastq.gz",
        zipped_fasta_reads = "data/reads/" + config["reads"] + ".fasta.gz",
        unzipped_fastq_reads = "data/reads/" + config["reads"] + ".fastq",
        unzipped_fasta_reads = "data/reads/" + config["reads"] + ".fasta"
    run:
        shell("""
        if [ ! -e {output[zipped_fastq_reads]} ] && [ ! -e {output[zipped_fasta_reads]} ] && [ ! -e {output[unzipped_fastq_reads]} ] && [ ! -e {output[unzipped_fasta_reads]} ]; then
            echo "Error: None of the expected input read files found. Rule failed. Please ensure your reads are in the data/reads/ directory."
            exit 1  # Exit with a non-zero status to indicate failure
        fi
        """)




#         # Define the Snakemake rule for downloading the NT database
# rule download_nt_db:
#     output:
#         os.path.join(config["ncbi_nt_path"], "nt.115.nin")
#     params:
#         ncbi_path = config["ncbi_nt_path"],
        # log = f"{config['assembly']}/logs/{rule}.log",
#     shell:
#         """
#         cd {params[ncbi_path]}
#         update_blastdb.pl --passive --decompress nt
#         cd - 2> {params[log]}
#         """
# blast database runs to nt.115. around 350 Gb of storage required

#  # Download the taxdb archive
# perl update_blastdb.pl taxdb
# # Install it in the BLASTDB directory
# gunzip -cd taxdb.tar.gz | (cd $BLASTDB; tar xvf - )

# Generate individual files for each scaffold
# rule splinter_assembly:
#     input:
#         assembly = "data/assemblies/" + config["assembly"] + ".fasta",
#         scaffolds = all_scaffs
#     output:
#         config["assembly"] + "/outputs/scaffolds/{all_scaffs}.fasta",
#     params:
#         scaffolds = config["assembly"] + "/outputs/scaffolds/",
        # log = f"{config['assembly']}/logs/{rule}.log",
#     shell:
#         """
#         cat {input} | awk '{{if (substr($0, 1, 1)=='>') {{filename=(substr($0,2) '.fasta'}} print $0 >> {params[scaffolds]}/filename
#         close({params[scaffolds]}/filename)}} 2> {params[log]}'
#         """

# # Unzip reads if zipped
# rule zip_fastq_to_fasta:
#     input:
#         reads = "data/reads/" + config["reads"] + ".fastq.gz",
#     output:
#         reads = "data/reads/" + config["reads"] + ".fasta",
#     params:
#         log = f"{config['assembly']}/logs/{rule}.log",
#     shell:
#         "zcat -c {input[reads]} | seqkit fq2fa | cat > {output} 2> {params[log]}"

# # Unzip reads if zipped
# rule zip_fasta_to_fasta:
#     input:
#         reads = "data/reads/" + config["reads"] + ".{ext}",
#     output:
#         reads = "data/reads/" + config["reads"] + ".fasta",
#     run:
#         if "fastq.gz" in input[0].path:
#             shell("echo 'Processing as FASTQ' > output.txt")
#             # Add your command for processing FASTQ files here
#         elif "fasta.gz" in input[0].path:
#             shell("echo 'Processing as FASTA' > output.txt")
#             # Add your command for processing FASTA files here
#         else:
#             raise ValueError("Unsupported reads file extension. Please use either fastq.gz or fasta.gz.")
