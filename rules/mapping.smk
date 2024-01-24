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
        seq_tech = "map-" + config["seq_tech"],
        log = f"{config['assembly']}/logs/mapping.log",
    shell:
        "minimap2 -t {params[threads]} -ax {params[seq_tech]} {input[assembly]} {input[reads]} > {output} 2> {params[log]}"  


# Convert the output sam file into a bam
rule conversion:
    input:
        sam = f"{config['assembly']}/outputs/mapping/{config['reads']}.sam"
    output:
        bam = f"{config['assembly']}/outputs/mapping/{config['reads']}.bam"
    params:
        log = f"{config['assembly']}/logs/conversion.log",
    shell:
        "samtools view -b -S {input[sam]} > {output[bam]} 2> {params[log]}"


# Sort the output bam file
rule sorting:
    input:
        bam = f"{config['assembly']}/outputs/mapping/{config['reads']}.bam"
    output:
        sorted_bam = f"{config['assembly']}/outputs/mapping/{config['reads']}.sorted.bam"
    params:
        log = f"{config['assembly']}/logs/sorting.log",
    shell:
        "samtools sort {input[bam]} > {output[sorted_bam]} 2> {params[log]}"


# Index the sorted bam file
rule samtools_index:
    input:
       bam = f"{config['assembly']}/outputs/mapping/{config['reads']}.sorted.bam"
    output:
       bai = f"{config['assembly']}/outputs/mapping/{config['reads']}.sorted.bam.bai",
    params:
        log = f"{config['assembly']}/logs/samtools_index.log",
    shell:
       "samtools index {input[bam]} > {output[bai]} 2> {params[log]}"
