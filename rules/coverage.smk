# COVERAGE RELATED THINGS

# Run mosdepth on the assembly 
rule mosdepth:
    input:
        bam = f"{config['assembly']}/outputs/mapping/{config['reads']}.sorted.bam",
        bai = f"{config['assembly']}/outputs/mapping/{config['reads']}.sorted.bam.bai",
    output:
        f"{config['assembly']}/reports/coverage/mosdepth/{config['assembly']}.mosdepth.summary.txt",
        f"{config['assembly']}/reports/coverage/mosdepth/{config['assembly']}.mosdepth.global.dist.txt"
    params:
        threads = config["threads"],
        out_pfx = f"{config['assembly']}/reports/coverage/mosdepth/{config['assembly']}"
    shell:
        "mosdepth --threads {params[threads]} {params[out_pfx]} {input[bam]}"


# Plot the output of mosdepth
rule mosdepth_plots:
    input:
        dist = f"{config['assembly']}/reports/coverage/mosdepth/{config['assembly']}.mosdepth.global.dist.txt"
    output:
        report(
            f"{config['assembly']}/reports/coverage/mosdepth/{config['assembly']}.dist.html",
            caption="../docs/captions/mosdepth.rst",
            category="Coverage analysis"
        )
    shell:
        "python3 scripts/plot_dist.py -o {output} {input[dist]}"


# Calculate base specific assembly coverage with samtools
rule get_coverage:
    input:
        bam = f"{config['assembly']}/outputs/mapping/{config['reads']}.sorted.bam",
        bai = f"{config['assembly']}/outputs/mapping/{config['reads']}.sorted.bam",
    output:
        f"{config['assembly']}/reports/coverage/{config['assembly']}.coverage"
    shell:
        "samtools depth {input[bam]} > {output}"


# # Break the coverage file into individual scaffold files
# rule scaffold_coverage:
#     input:
#         assembly = f"data/assemblies/{config['assembly']}.fasta"
#     output:
#         f"{config['assembly']}/reports/coverage/{all_scaffs}.coverage"
#     shell:
#         "awk '$1 == {all_scaffs} '{{print $0}}' {input} > {output}"


# Generate plots of coverage across the assembly

rule assembly_coverage_plot:
    input:
        coverage = f"{config['assembly']}/reports/coverage/{config['assembly']}.coverage"
    output:
        f"{config['assembly']}/reports/coverage/{config['assembly']}.coverage.png"
    script:
        "../scripts/assembly_coverage.R"


# # Generate plots of coverage across the scaffolds
# rule scaffold_coverage_plots:
#     input:
#         f"{config['assembly']}/reports/coverage/{all_scaffs}.coverage"
#     output:
#         f"{config['assembly']}/reports/coverage/{all_scaffs}.coverage.png"
#     script:
#         "../scripts/scaffold_coverage.R"
