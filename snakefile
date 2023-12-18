# assembly appraisal 

report: "docs/captions/workflow.rst"

configfile: "config.yaml"

include: "rules/coverage.smk"
#include: "rules/blobplots.smk"
# include: "rules/merqury.smk"
#include: "rules/characterisation.smk"
#include: "rules/completeness.smk"
#include: "rules/mito.smk"
# include: "rules/checks_and_transformations.smk"
include: "rules/mapping.smk"
#include: "rules/nucmer.smk"
#include: "rules/quast.smk"
#include: "rules/pair_analysis.smk"
#include: "rules/variant_calling.smk"

###############################################################################

# TODO fix this rule below me to get the depth plots working
# rule get_scaffs:
#     output:
#         all_scaffs
#     run:
#         all_scaffs = []
#         with open(snakemake.input[assembly], "r") as f:
#                 for record in SeqIO.parse(f, "fasta"):
#                         all_scaffs.append(record.description)


#########################################



rule final_outputs:
    input:
# GENOME PROFILING
#        genomescope = config["assembly"] + "/reports/genomescope/plot.png",
        # smudgeplot = config["assembly"] + "/reports/smudge/smudgeplot_smudgeplot.png",
# BLOBPLOTS
#        blob_table = config["assembly"] + "/reports/blobtools/" + config["assembly"] + ".blobDB.table.txt",
#        blob_plot = config["assembly"] + "/reports/blobtools/" + config["assembly"] + ".blobDB.json.bestsum.phylum.p8.span.100.blobplot.bam0.png",
# NUCMER
#        self_v_self = config["assembly"] + "/outputs/nucmer/nucmer.self_v_self.png",
#        self_v_ref = config["assembly"] + "/outputs/nucmer/nucmer.self_v_ref.png",
#         # dna_diff =
# PAIRS ANALYSIS
        # blast_pairs = config["assembly"] + "/reports/pairs_analysis/blast/blast.onlyPairs.tsv",
        # dotplots = directory(config["assembly"] + "/reports/pairs_analysis/nucmer/pairs"),
#         # dna_diff =
# QUAST
#        quast_report = config["assembly"] + "/reports/quast/report.html",
# CEGMA
#        completeness_report = config["assembly"] + "/reports/cegma/" + config["assembly"] + ".completeness_report",
# COVERAGE
        plots_initial = config["assembly"] + "/reports/coverage/mosdepth/initial_" + config["assembly"] + ".dist.html",
# MITO
#        mito_tagged = config["assembly"] + "/outputs/assemblies/" + config["assembly"] + ".mito_tagged.fasta",
#        no_mito = config["assembly"] + "/outputs/assemblies/" + config["assembly"] + ".no_mito.fasta",
# VARIANT CALLING
#        sniffles = config["assembly"] + "/outputs/variant_calling/" + config["assembly"] + ".vcf",

# MERQURY 
#         # merqury_mrls =
#         # merqury_out =
# LOGS
#        config = config["assembly"] + "/logs/config.log",
#        environment = config["assembly"] + "/logs/environment.log"



rule log_config:
    input:
        "config.yaml"
    output:
        config = config["assembly"] + "/logs/config.log",
    shell:
        "cp {input} {output}"

rule log_environment:
    output:
        config = config["assembly"] + "/logs/environment.log",
    shell:
        "conda env export -f {output}"


