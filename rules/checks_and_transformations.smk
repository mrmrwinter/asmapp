# INPUT CHECKS AND INITIAL DATA TRANSFORMATIONS

# Check for presence of input assembly
rule input_assembly:
    output:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",

# Unzip reads if zipped
rule zip_fastq_to_fasta:
    input:
        reads = "data/reads/" + config["reads"] + ".fastq.gz",
    output:
        reads = "data/reads/" + config["reads"] + ".fasta",
    shell:
        "zcat -c {input[reads]} | seqkit fq2fa | cat > {output}"


