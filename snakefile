# assembly appraisal 

report: "workflow.rst"

configfile: "config.yaml"

include: "rules/pair_alignment.smk"
include: "rules/coverage.smk"
include: "rules/blobplots.smk"
# include: "rules/merqury.smk"
include: "rules/characterisation.smk"
include: "rules/cegma.smk"
include: "rules/mito.smk"
include: "rules/checks_and_transformations.smk"
include: "rules/mapping.smk"
include: "rules/nucmer.smk"
# include: "rules/karyon.smk"
include: "rules/quast.smk"
include: "rules/pair_analysis.smk"


include: "rules/variant_calling.smk"

###############################################################################

rule all:
    input:
# INPUT CHECK
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",
        reads = "data/reads/" + config["reads"] + ".fastq.gz",
# TAGGING
        # intial_tagged = config["assembly"] + "/outputs/cegma/tagged_initial_assembly.fasta",
# GENOME PROFILING
        # genomescope # TODO
        smudgeplot = config["assembly"] + "/reports/smudge/smudgeplot_smudgeplot.png",
# BLOBPLOTS
        blob_table = config["assembly"] + "/reports/blobtools/" + config["assembly"] + ".blobDB.table.txt",
        blob_plot = config["assembly"] + "/reports/blobtools/" + config["assembly"] + ".blobDB.json.bestsum.phylum.p8.span.100.blobplot.bam0.png",
# NUCMER
        self_v_self = config["assembly"] + "/outputs/nucmer/nucmer.self_v_self.png",
        self_v_ref = config["assembly"] + "/outputs/nucmer/nucmer.self_v_ref.png",
#         # dna_diff =
# # # BLAST TABLES
# #         tsv = config["assembly"] + "/outputs/blast/blast.out",
# #         pairs_tsv = config["assembly"] + "/outputs/blast/blast.out",
        blast_pairs = config["assembly"] + "/outputs/blast/blast.onlyPairs.tsv",
# # PAIRS ALIGNMENT
#         # dotplots =
# #         # dna_diff =
# # # QUAST
# #         # busco_lib =
        quast_report = config["assembly"] + "/reports/quast/report.txt",
        # TODO replace QUAST with my own installs and scripts
# # MERQURY # TODO add merqury
#         # merqury_mrls =
#         # merqury_out =
# # CEGMA
        completeness_report = config["assembly"] + "/reports/cegma/" + config["assembly"] + ".completeness_report",
# TODO add BUSCO v5
# # COVERAGE
        plots_initial = config["assembly"] + "/reports/coverage/mosdepth/initial_" + config["assembly"] + ".dist.html",
# # MITO
        mito_tagged = config["assembly"] + "/outputs/assemblies/" + config["assembly"] + ".mito_tagged.fasta",
        # no_mito = config["assembly"] + "/outputs/assemblies/" + config["assembly"] + ".no_mito.fasta",
# # VARIANT CALLING
        sniffles = config["assembly"] + "outputs/variant_calling/" + config["assembly"] + ".vcf",






