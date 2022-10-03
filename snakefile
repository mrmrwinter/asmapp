# ASMAPP - assembly appraisal workflow

report: "reports/workflow.rst"

configfile: "config.yaml"

include: [
        "rules/coverage.smk", 
        "rules/blobplots.smk",
        # "rules/merqury.smk",
        "rules/characterisation.smk",
        "rules/cegma.smk",
        "rules/mito.smk",
        "rules/checks_and_transformations.smk",
        "rules/mapping.smk",
        "rules/nucmer.smk",
        "rules/quast.smk",
        "rules/pair_analysis.smk",
        "rules/variant_calling.smk"
]

###############################################################################

rule all:
    input:
# GENOME PROFILING
        smudgeplot = config["assembly"] + "/reports/smudge/smudgeplot_smudgeplot.png",
# BLOBPLOTS
        blob_table = config["assembly"] + "/reports/blobtools/" + config["assembly"] + ".blobDB.table.txt",
        blob_plot = config["assembly"] + "/reports/blobtools/" + config["assembly"] + ".blobDB.json.bestsum.phylum.p8.span.100.blobplot.bam0.png",
# NUCMER
        self_v_self = config["assembly"] + "/outputs/nucmer/nucmer.self_v_self.png",
        self_v_ref = config["assembly"] + "/outputs/nucmer/nucmer.self_v_ref.png",
#         # dna_diff =
# PAIRS ANALYSIS
        blast_pairs = config["assembly"] + "/outputs/blast/blast.onlyPairs.tsv",
        dotplots = directory(config["assembly"] + "reports/nucmer/pairs"),
#         # dna_diff =
# QUAST
        quast_report = config["assembly"] + "/reports/quast/report.html",
# CEGMA
        completeness_report = config["assembly"] + "/reports/cegma/" + config["assembly"] + ".completeness_report",
# COVERAGE
        plots_initial = config["assembly"] + "/reports/coverage/mosdepth/initial_" + config["assembly"] + ".dist.html",
# MITO
        mito_tagged = config["assembly"] + "/outputs/assemblies/" + config["assembly"] + ".mito_tagged.fasta",
        no_mito = config["assembly"] + "/outputs/assemblies/" + config["assembly"] + ".no_mito.fasta",
# VARIANT CALLING
        sniffles = config["assembly"] + "/outputs/variant_calling/" + config["assembly"] + ".vcf",

# MERQURY 
#         # merqury_mrls =
#         # merqury_out =

# # TODO add BUSCO v5
#         # TODO replace QUAST with my own installs and scripts
# TODO add merqury
# genomescope # TODO