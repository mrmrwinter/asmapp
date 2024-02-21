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
        out_pfx = f"{config['assembly']}/reports/coverage/mosdepth/{config['assembly']}",
        log = f"{config['assembly']}/logs/mosdepth.log",
    shell:
        "mosdepth --threads {params[threads]} {params[out_pfx]} {input[bam]} 2> {params[log]}"


# Plot the output of mosdepth
rule mosdepth_plots:
    input:
        dist = f"{config['assembly']}/reports/coverage/mosdepth/{config['assembly']}.mosdepth.global.dist.txt"
    output:
        # html = f"{config['assembly']}/reports/coverage/mosdepth/{config['assembly']}.dist.html",
        report = report(
            f"{config['assembly']}/reports/coverage/mosdepth/{config['assembly']}.dist.html",
            caption="../docs/captions/mosdepth.rst",
            category="Coverage analysis"
        )
    params:
        log = f"{config['assembly']}/logs/mosdepth_plots.log",
    shell:
        "python3 scripts/plot_dist.py -o {output} {input[dist]} 2> {params[log]}"


# Calculate base specific assembly coverage with samtools
rule get_coverage:
    input:
        bam = f"{config['assembly']}/outputs/mapping/{config['reads']}.sorted.bam",
        bai = f"{config['assembly']}/outputs/mapping/{config['reads']}.sorted.bam",
    output:
        f"{config['assembly']}/reports/coverage/{config['assembly']}.coverage"
    params:
        log = f"{config['assembly']}/logs/get_coverage.log",
    shell:
        "samtools depth {input[bam]} > {output} 2> {params[log]}"


# # Generate individual files for each scaffold
# all_scaffs = []

# with open(f"{config['assembly']}/logs/scaffold_list.txt", "r") as file:
#     # Read each line from the file
#     for line in file:
#         # Strip any leading/trailing whitespace and append the line to the list
#         all_scaffs.append(line.strip())


# # Break the coverage file into individual scaffold files
# rule scaffold_coverage:
#     input:
#         assembly_coverage = f"{config['assembly']}/reports/coverage/{config['assembly']}.coverage",
#         scaffolds = f"{config['assembly']}/logs/scaffold_list.txt"
#     output:
#         # directory(config['assembly'] + "/reports/coverage/scaffolds")
#         expand(config['assembly'] + "/reports/coverage/scaffolds/{scaffold}.coverage", scaffold = all_scaffs)
#     params:
#         out_pfx = config['assembly'] + "/reports/coverage/scaffolds"
#     shell:
#         """
#         mkdir -p {params.out_pfx}

#         awk '{{print > "{params.out_pfx}/"$1".coverage"}}' {input.assembly_coverage}
#         """

# rule scaffold_coverage_plot:
#     input:
#         expand(config['assembly'] + "/reports/coverage/scaffolds/{scaffold}.coverage", scaffold = all_scaffs)
#     output:
#         files = expand(config['assembly'] + "/reports/coverage/scaffolds/plots/{scaffold}.coverage.png", scaffold = all_scaffs),
#         path = directory(f"{config['assembly']}/reports/coverage/scaffolds/plots")
#     script:
#         "../scripts/coverage_plots.py"

