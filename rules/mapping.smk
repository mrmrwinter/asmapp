###############################################################################

# MAPPING AND TRANSFORMATIONS FOR initial ASSEMBLY

rule initial_mapping:
    input:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",
        reads = "data/reads/" + config["reads"] + ".fastq.gz"
    output:
        sam = config["assembly"] + "/outputs/initial/initial_asm.sam"
    params:
        threads = config["threads"],
        seq_tech = "map-" + config["seq_tech"]
    shell:
        "minimap2 -t {params[threads]} -ax {params[seq_tech]} {input[assembly]} {input[reads]} > {output}"


rule initial_conversion:
    input:
        sam = config["assembly"] + "/outputs/initial/initial_asm.sam"
    output:
        bam = config["assembly"] + "/outputs/initial/initial_asm.bam"
    shell:
        "samtools view -b -S {input} > {output}"


rule initial_sorting:
    input:
        config["assembly"] + "/outputs/initial/initial_asm.bam"
    output:
        bam = config["assembly"] + "/outputs/initial/initial_asm.sorted.bam"
    shell:
        "samtools sort {input} > {output}"


rule initial_samtools_index:
    input:
       bam = config["assembly"] + "/outputs/initial/initial_asm.sorted.bam"
    output:
       bai = config["assembly"] + "/outputs/initial/initial_asm.sorted.bam.bai",
    shell:
       "samtools index {input[bam]} > {output[0]}"



###############################################################################