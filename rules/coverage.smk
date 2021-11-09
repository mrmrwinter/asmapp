# COVERAGE RELATED THINGS




rule mosdepth:
    input:
        initial = "data/assemblies/" + config["assembly"] + ".fasta",
        reads = "data/reads/" + config["reads"] + ".fastq.gz"
    output:
        config["assembly"] + "reports/coverage/mosdepth/" + config["assembly"] + "mosdepth.summary.txt"
    params:
        threads = config["threads"],
        out_pfx = config["assembly"] + "reports/coverage/mosdepth/" + config["assembly"]
    shell:
        "mosdepth --threads {params[threads]} {params[out_pfx]}"
