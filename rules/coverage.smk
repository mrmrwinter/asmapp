# COVERAGE RELATED THINGS

# Run mosdepth on the assembly 
rule mosdepth:
    # conda:
    #     "../envs/coverage.yaml"
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

# Plot the output of mosdepth
rule mosdepth_plots:
    # conda:
    #     "../envs/coverage.yaml"
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

# Calculate base specific assembly coverage with samtools
rule get_coverage:
    input:
        bam = config["assembly"] + "/outputs/initial/initial_asm.sorted.bam",
        bai = config["assembly"] + "/outputs/initial/initial_asm.sorted.bam.bai",
    output:
        config["assembly"] + "/reports/coverage/" + config["assembly"] + ".coverage"
    shell:
        "samtools depth {input[bam]} > {output}"

# Break the coverage file into individual scaffold files
rule scaffold_coverage:
    input:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta"
    output:
        config["assembly"] + "/reports/coverage/{all_scaffs}.coverage"
    shell:
        "awk '$1 == {all_scaffs} '{{print $0}}' {input} > {output}"

# Generate plots of coverage across the assembly
rule assembly_coverage_plot:
    input:
        config["assembly"] + "/reports/coverage/" + config["assembly"] + ".coverage"
    output:
        config["assembly"] + "/reports/coverage/" + config["assembly"] + ".coverage.png"
    script:
        "../scripts/assembly_coverage.R"

# Generate plots of coverage across the scaffolds
rule scaffold_coverage_plots:
    input:
        config["assembly"] + "/reports/coverage/{all_scaffs}.coverage"
    output:
        config["assembly"] + "/reports/coverage/{all_scaffs}.coverage.png"
    script:
        "../scripts/scaffold_coverage.R"
