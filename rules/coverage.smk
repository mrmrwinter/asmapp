# COVERAGE RELATED THINGS

rule mosdepth:
    conda:
        "../envs/coverage.yaml"
    input:
        bam = config["assembly"] + "/outputs/initial/initial_asm.sorted.bam",
        bai = config["assembly"] + "/outputs/initial/initial_asm.sorted.bam.bai",
    output:
        config["assembly"] + "/reports/coverage/mosdepth/initial_" + config["assembly"] + ".mosdepth.summary.txt",
        config["assembly"] + "/reports/coverage/mosdepth/initial_" + config["assembly"] + ".mosdepth.global.dist.txt"
    params:
        threads = config["threads"],
        out_pfx = config["assembly"] + "/reports/coverage/mosdepth/initial_" + config["assembly"]
    shell:
        "mosdepth --threads {params[threads]} {params[out_pfx]} {input[0]}"

rule mosdepth_plots:
    conda:
        "../envs/coverage.yaml"
    input:
        config["assembly"] + "/reports/coverage/mosdepth/initial_" + config["assembly"] + ".mosdepth.global.dist.txt"
    output:
        report(
            config["assembly"] + "/reports/coverage/mosdepth/initial_" + config["assembly"] + ".dist.html",
            caption="../docs/captions/mosdepth.rst",
            category="Coverage analysis"
        )
    shell:
        "python3 scripts/plot_dist.py -o {output} {input}"
