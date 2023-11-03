# INPUT CHECKS AND INITIAL DATA TRANSFORMATIONS

# Check for presence of input assembly
rule input_assembly:
    output:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",

# Check for presenc of input reads
rule input_reads:
    output:
        reads = "data/reads/" + config["reads"] + ".fastq.gz",




#         # Define the Snakemake rule for downloading the NT database
rule download_nt_db:
    output:
        os.path.join(config["ncbi_nt_path"], "nt.115.nin")
    params:
        config["ncbi_nt_path"]
    shell:
        """
        cd {params}
        update_blastdb.pl --passive --decompress nt
        cd -
        """
# blast database runs to nt.115. around 350 Gb of storage required


# Generate individual files for each scaffold
# rule splinter_assembly:
#     input:
#         assembly = "data/assemblies/" + config["assembly"] + ".fasta",
#         scaffolds = all_scaffs
#     output:
#         config["assembly"] + "/outputs/scaffolds/{all_scaffs}.fasta",
#     params:
#         config["assembly"] + "/outputs/scaffolds/"
#     shell:
#         """
#         cat {input} | awk '{{if (substr($0, 1, 1)=='>') {{filename=(substr($0,2) '.fasta'}} print $0 >> {params}/filename
#         close({params}/filename)}}'
#         """

# Unzip reads if zipped
rule reads_to_fasta:
    input:
        reads = "data/reads/" + config["reads"] + ".fastq.gz",
    output:
        reads = "data/reads/" + config["reads"] + ".fasta",
    shell:
        "zcat -c {input} | seqkit fq2fa | cat > {output}"
