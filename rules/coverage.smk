# COVERAGE RELATED THINGS

rule mosdepth:
    input:
        initial = "data/assemblies/" + config["assembly"] + ".fasta",
        reads = "data/reads/" + config["reads"] + ".fastq.gz"
    output:
        
    shell:
        "mosdepth"
