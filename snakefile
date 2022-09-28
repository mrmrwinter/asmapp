# assembly appraisal (karyinon workflow)

report: "workflow.rst"

configfile: "config.yaml"

include: "rules/pair_alignment.smk"
include: "rules/coverage.smk"
include: "rules/blobplots.smk"
# include: "rules/merqury"
include: "rules/characterisation.smk"
include: "rules/cegma.smk"
# include: "rules/mito.smk"
include: "rules/checks_and_transformations.smk"
# include: "rules/mapping.smk"
# include: "rules/nucmer.smk"
# include: "rules/karyon.smk"

###############################################################################

rule all:
    input:
# INPUT CHECK
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",
        reads = "data/reads/" + config["reads"] + ".fastq.gz",
# TAGGING
        # intial_tagged = config["assembly"] + "/outputs/cegma/tagged_initial_assembly.fasta",
        # collapsed_tagged = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta",
# GENOME PROFILING
#        smudgeplot = config["assembly"] + "/reports/smudge/smudgeplot_smudgeplot.png",
# BLOBPLOTS
        blob_table = config["assembly"] + "/reports/blobtools/" + config["assembly"] + ".blobDB.table.txt",
        blob_plot = config["assembly"] + "/reports/blobtools/" + config["assembly"] + ".blobDB.json.bestsum.phylum.p8.span.100.blobplot.bam0.png",
# # WHOLE GENOME ALIGNMENTS
        nucmer = config["assembly"] + "/reports/nucmer/nucmer.initial.delta",
        initial = config["assembly"] + "/reports/nucmer/nucmer.initial.png",
        # initial_v_collapsed = config["assembly"] + "/reports/nucmer/nucmer.initial_v_collapsed.png",
        # collapsed_v_ref = config["assembly"] + "/reports/nucmer/nucmer.collapsed_v_ref.delta",
        initial_v_ref = config["assembly"] + "/reports/nucmer/nucmer.initial_v_ref.delta",
#         # dna_diff =
# # # BLAST TABLES
# #         tsv = config["assembly"] + "/reports/blast/blast.out",
#         only_pairs = config["assembly"] + "/reports/blast/blast.onlyPairs.tsv",
# #         initial_tsv = config["assembly"] + "/reports/blast/initial_blast.out",
#         initial_only_pairs = config["assembly"] + "/reports/blast/initial_blast.onlyPairs.tsv",
# # PAIRS ALIGNMENT
#         # dotplots =
# #         # dna_diff =
# # # QUAST
# #         # busco_lib =
        quast_report = config["assembly"] + "/reports/quast/report.txt",
# # MERQURY
#         # merqury_mrls =
#         # merqury_out =
# # CEGMA
        completeness_report = config["assembly"] + "/outputs/cegma/" + config["assembly"] + ".completeness_report",
# # COVERAGE
        # mosdepth_haplome = config["assembly"] + "/reports/coverage/mosdepth/collapsed_" + config["assembly"] + ".mosdepth.summary.txt",
        # plots_haplome = config["assembly"] + "/reports/coverage/mosdepth/collapsed_" + config["assembly"] + ".dist.html",
        plots_initial = config["assembly"] + "/reports/coverage/mosdepth/initial_" + config["assembly"] + ".dist.html",
# # KARYON
#         # flagstats = config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.flagstat",
#         vcf = config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.vcf",
#         # mpileup = config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.mpileup",
#         plot = config["assembly"] + "/outputs/plots/plot.png",
# # MITO
        # mito_tagged = "data/assemblies/" + config["assembly"] + ".mito_tagged.fasta",
        # no_mito = "data/assemblies/" + config["assembly"] + ".no_mito.fasta",








###############################################################################
# FLUFF
##################################################

# rule dotplots:
#     input:
#         tsv = config["assembly"] + "/reports/mapped_nonself_hits.sam",
#         assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta"
#     output:
#         config["assembly"] + "/reports/dotplots/dotplot.png"
#     params:
#         out_dir = config["assembly"] + "/reports/dotplots/"
#     run:
#         from Bio import SeqIO
#         import pandas as pd
#         import os
#         import sys
#         from readpaf import parse_paf
#         import dotplot
#
#         shell("mkdir -p tmp/")
#
#         tmp_tbl = input[0].replace(".paf",".cleaned.sam")
#         print(tmp_tbl)
#         shell("sed '/^@/ d' < " + input[0] + " > " + tmp_tbl)
#
#         df = pd.read_csv(tmp_tbl, sep = '\t')
#         seqs = input[1]
# #
#         for index, value in df.iterrows():
#             query = value[0]
#             hit = value[2]
#             print(query + " pairs with " + hit)
#             for record in SeqIO.parse(seqs,'fasta'):
#                 q_seq=[]
#                 h_seq=[]
#                 if record.id == query:
#                     q_seq = ">" + record.id + "\n" + record.seq + "\n"
#                     with open('tmp/' + query + '.fasta', 'a+', newline='\n') as g:
#                         SeqIO.write(record, g, 'fasta-2line')
#                     g.close()
#                     print(q_seq)
#                 if record.id == hit:
#                     h_seq = ">" + record.id + "\n" + record.seq + "\n"
#                     with open('tmp/' + hit + '.fasta', 'a+', newline='\n') as r:
#                         SeqIO.write(record, r, 'fasta-2line')
#                     r.close()
#                     print(h_seq)
#             shell("dotplot --drawer matplotlib --fasta tmp/" + str(query) + ".fasta tmp/" + str(hit) + ".fasta > " + params[0] + str(query) + ".dotplot.png")



# run:
#     from Bio import SeqIO
#     import os
#
#     if sys.version_info[0] < 3:
#         from StringIO import StringIO
#     else:
#         from io import StringIO
#
#     for seq_record in SeqIO.parse(input[0], "fasta"):
#         contig = '>' + seq_record.description + '\n' + str(seq_record.seq)
# -negative_seqidlist exclude_me
# "blastn -query {input[0]} -db {params[1]} -outfmt '6 qseqid sseqid pident' -out {params[0]}/blast.out"
#                 # q =
#                 # h = df.[1]
#                 # def get_first_pd(condition, df):
#                 #     return df[condition(df)].iloc[0]
#                 # q = str(df.iloc[1]).split('\\t')[0]
#                 # q = q.replace('redundansscaffold','')
#                 # q = int(q)
#                 # h = str(df.iloc[1]).split('\\t')[0]
#                 # h = h.replace('redundansscaffold','')
#                 # h = int(h)
#                 # # hit_row = get_first_pd(lambda x: x(contig != h), df)
#                 # print(h)
#                 # writer.writerow(hit_row)
#
#                 shell("rm " + seq_record.id + "*.out")
#
#         g.close()

#
#

# rule mapping_back:
#     input:
#         assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta"
#     output:
#         tsv = config["assembly"] + "/reports/minimap2/mapped_nonself_hits.sam"
#     params:
#         threads = config["threads"]
#     shell:
#         "minimap2 -ax asm5 -D {input[0]} {input[0]} > {output}"
#
# # rule removing_nonself:
# #     input:
# #         tsv = config["assembly"] + "/reports/mapped_hits.sam"
# #     output:
# #         tsv = config["assembly"] + "/reports/mapped_nonself_hits.sam"
# #     shell:
# #         "samtools view -F0x900 {input} > {output}"
