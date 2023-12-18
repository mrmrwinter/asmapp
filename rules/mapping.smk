# MAPPING AND TRANSFORMATIONS 


# Map the reads to the assembly
rule mapping:
    input:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",
        reads = "data/reads/" + config["reads"] + ".fastq.gz"
    output:
        sam = f"{config['assembly']}/outputs/mapping/{config['reads']}.sam"
    params:
        threads = config["threads"],
        seq_tech = "map-" + config["seq_tech"]
    shell:
        "minimap2 -t {params[threads]} -ax {params[seq_tech]} {input[assembly]} {input[reads]} > {output}"  

# Convert the output sam file into a bam
rule conversion:
    input:
        sam = f"{config['assembly']}/outputs/mapping/{config['reads']}.sam"
    output:
        bam = f"{config['assembly']}/outputs/mapping/{config['reads']}.bam"
    shell:
        "samtools view -b -S {input} > {output}"

# Sort the output bam file
rule sorting:
    input:
        f"{config['assembly']}/outputs/mapping/{config['reads']}.bam"
    output:
        bam = f"{config['assembly']}/outputs/mapping/{config['reads']}.sorted.bam"
    shell:
        "samtools sort {input} > {output}"

# Index the sorted bam file
rule samtools_index:
    input:
       bam = f"{config['assembly']}/outputs/mapping/{config['reads']}.sorted.bam"
    output:
       bai = f"{config['assembly']}/outputs/mapping/{config['reads']}.sorted.bam.bai",
    shell:
       "samtools index {input[bam]} > {output[0]}"
