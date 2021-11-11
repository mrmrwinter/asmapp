# COVERAGE RELATED THINGS


rule mosdepth_collapsed:
    conda:
        "../envs/coverage.yaml"
    input:
        bam = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sorted.bam",
        bai = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sorted.bam.bai",
    output:
        config["assembly"] + "/reports/coverage/mosdepth/collapsed_" + config["assembly"] + ".mosdepth.summary.txt",
        config["assembly"] + "/reports/coverage/mosdepth/collapsed_" + config["assembly"] + ".mosdepth.global.dist.txt"
    params:
        threads = config["threads"],
        out_pfx = config["assembly"] + "/reports/coverage/mosdepth/collapsed_" + config["assembly"]
    shell:
        "mosdepth --threads {params[threads]} {params[out_pfx]} {input[0]}"

rule mosdepth_collapsed_plots:
    conda:
        "../envs/coverage.yaml"
    input:
        config["assembly"] + "/reports/coverage/mosdepth/collapsed_" + config["assembly"] + ".mosdepth.global.dist.txt"
    output:
        config["assembly"] + "/reports/coverage/mosdepth/collapsed_" + config["assembly"] + ".dist.html"
    shell:
        "python3 scripts/plot_dist.py -o {output} {input}"
