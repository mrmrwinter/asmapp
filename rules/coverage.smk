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


rule get_coverage:
    input:
        bam = config["assembly"] + "/outputs/initial/initial_asm.sorted.bam",
        bai = config["assembly"] + "/outputs/initial/initial_asm.sorted.bam.bai",
    output:
        config["assembly"] + "/reports/coverage/" + config["assembly"] + ".coverage"
    shell:
        "samtools depth {input[bam]} > {output}"


# TODO get the following rule working
rule scaffold_coverage:
    input:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta"
    output:
        config["assembly"] + "/reports/coverage/{all_scaffs}.coverage"
    shell:
        "awk '$1 == {all_scaffs} '{{print $0}}' {input} > {output}"


rule assembly_coverage_plot:
    input:
        config["assembly"] + "/reports/coverage/" + config["assembly"] + ".coverage"
    output:
        config["assembly"] + "/reports/coverage/" + config["assembly"] + ".coverage.png"
    script:
        "../scripts/assembly_coverage.R"


rule scaffold_coverage_plots:
    input:
        config["assembly"] + "/reports/coverage/{all_scaffs}.coverage"
    output:
        config["assembly"] + "/reports/coverage/{all_scaffs}.coverage.png"
    script:
        "../scripts/scaffold_coverage.R"