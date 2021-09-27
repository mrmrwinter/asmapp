# PAIRS ALIGNMENTS FROM NUCMER

# rule pair_alignment:
#     input:
#        blast =
#     output:
#         only_pairs_table =
#     run:
#         import pandas as pd
#
#         blast_output = pd.read_csv(snakemake.input[0], sep="\t", header = None) # snakemake.input[0] is the blast table
#
#         initial_pairs = pd.DataFrame(columns = ['query', 'hit'])
#
#         for index, value in blast_output.iterrows():
#             if value[0] != value[1]:
#                 initial_pairs.loc[index, ['query']] = value[0]
#                 initial_pairs.loc[index, ['hit']] = value[1]
#
#         only_initial_pairs = initial_pairs.drop_duplicates()
#
#         only_initial_pairs.to_csv(snakemake.output[0], sep='\t')
#


# IS THERE A WAY TO TAKE A SAMPLE LIST FROM THE ONLY PAIRS COLUMN HERE?
# IF SO I WONT HAVE TO USE THE DIRECTORY OUTPUT FLAG WHICH IS A pair_alignment
# CAN INSTEAD LOOP THROUGH THE LIST FOR EACH INDIVIDUAL ASSEMBLY AS THE LIST OF PAIRS WILL BE NOVEL AND INDIVIDUAL

#
# rule single_scaffold_extraction:
#     input:
#         assembly =
#     output:
#         directory(config["assembly"] + "tmp/")
#     params:
#         out_dir = config["assembly"] + "tmp/"
#     shell:
#         """
#         mkdir -p {params[0]}
#         awk '/^>/ {OUT='{params[0]}' substr($0,2) ''.fasta'}; OUT{print >OUT}' {input[assembly]}
#         """
#
#
# # nucmer analysis/alignment
# rule :
#     input:
#         directory(config["assembly"] + "tmp/")
#         only_pairs_table =
#     output:
#
#     params:
#         out_dir = config["assembly"] + "tmp/"
#     run:
#         import glob
#         import os
#
#         # os.system("mkdir nucmer_initial_purged/")
#
#         pairs = pd.read_csv()
#
#         for index, value in only_initial_pairs.iterrows():
#             q = "tmp_initial_purged/" + value[0] + ".fasta"
#             h = "tmp_initial_purged/" + value[1] + ".fasta"
#             nucmer = "nucmer -p nucmer_initial_purged/nucmer." + str(index) + " " + q + " " + h
#             os.system(nucmer)
#
#
#         for delta in glob.glob("nucmer_initial_purged/*.delta"):
#             pair = delta.replace("nucmer_initial_purged/","").replace(".delta","")
#             mummer = "mummerplot -l -f --png --large " + delta + " -p nucmer_initial_purged/" + pair
#             os.system(mummer)
