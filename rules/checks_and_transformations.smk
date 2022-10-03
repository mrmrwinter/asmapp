# INPUT CHECKS

rule input_assembly:
    output:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",


rule input_reads:
    output:
        reads = "data/reads/" + config["reads"] + ".fastq.gz",


# # RENAME INITIAL CONTIGS
# rule initial_tagging:
#     input:
#     output:
#         assembly = "M_javanica_062022.final.renamed/outputs/cegma/tagged_initial_assembly.fasta"
#     run:
#         from Bio import SeqIO
#
#         to_add = "scaffold_"
#         with open(output[0], "w") as outputs:
#             for r in SeqIO.parse(input[0], "fasta"):
#                 r.id = (to_add + r.description).replace(" ", "_")
#                 r.description = r.id
#                 SeqIO.write(r, outputs, "fasta")
