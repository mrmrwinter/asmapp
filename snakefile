# assembly appraisal 

report: "docs/captions/workflow.rst"

configfile: "config.yaml"

include: "rules/coverage.smk"
include: "rules/blobplots.smk"
include: "rules/characterisation.smk"
include: "rules/completeness.smk"
include: "rules/mito.smk"
include: "rules/checks_and_transformations.smk"
include: "rules/mapping.smk"
include: "rules/nucmer.smk"
include: "rules/quast.smk"
#include: "rules/pair_analysis.smk"
include: "rules/variant_calling.smk"


rule final_outputs:
    input:
# GENOME PROFILING
        # genomescope = f"{config['assembly']}/reports/genomescope/plot.png",
        # smudgeplot = f"{config['assembly']}/reports/smudge/smudgeplot_smudgeplot.png",  # muted. requires big memory
# BLOBPLOTS
        # blob_table = f"{config['assembly']}/reports/blobtools/{config['assembly']}.blobDB.table.txt",
        # blob_plot = f"{config['assembly']}/reports/blobtools/{config['assembly']}.blobDB.json.bestsum.phylum.p8.span.100.blobplot.bam0.png",
# NUCMER
        self_v_self = f"{config['assembly']}/outputs/nucmer/nucmer.self_v_self.png",
        self_v_ref = f"{config['assembly']}/outputs/nucmer/nucmer.self_v_ref.png",
# PAIRS ANALYSIS
        # blast_pairs = f"{config['assembly']}/reports/pairs_analysis/blast/blast.onlyPairs.tsv",
        # dotplots = directory(f"{config['assembly']}/reports/pairs_analysis/nucmer/pairs"),
# QUAST
        quast_report = f"{config['assembly']}/reports/quast/report.html",
# CEGMA
        # completeness_report = f"{config['assembly']}/reports/cegma/{config['assembly']}.completeness_report",
        busco = f"{config['assembly']}/reports/busco/run_{config['busco_lineage']}/short_summary.txt",
# COVERAGE
        mosdepth_plot = f"{config['assembly']}/reports/coverage/mosdepth/{config['assembly']}.dist.html",
        all_coverage = f"{config['assembly']}/reports/coverage/{config['assembly']}.coverage",
# MITO
        mito_tagged = f"{config['assembly']}/outputs/assemblies/{config['assembly']}.mito_tagged.fasta",
        no_mito = f"{config['assembly']}/outputs/assemblies/{config['assembly']}.no_mito.fasta",
# VARIANT CALLING
        sniffles = f"{config['assembly']}/outputs/variant_calling/{config['assembly']}_{config['reads']}.vcf",
# LOGS
        config = f"{config['assembly']}/logs/config.log",
        environment = f"{config['assembly']}/logs/environment.log"



rule log_config:
    input:
        "config.yaml"
    output:
        config = f"{config['assembly']}/logs/config.log",
    shell:
        "cp {input} {output}"

rule log_environment:
    output:
        config = f"{config['assembly']}/logs/environment.log",
    shell:
        "conda env export -f {output}"


