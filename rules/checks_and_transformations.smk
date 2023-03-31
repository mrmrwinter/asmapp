# INPUT CHECKS


rule input_assembly:
    output:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",


rule input_reads:
    output:
        reads = "data/reads/" + config["reads"] + ".fastq.gz",


rule splinter_assembly:
    input:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",
        scaffolds = all_scaffs
    output:
        config["assembly"] + "/outputs/scaffolds/{all_scaffs}.fasta",
    params:
        config["assembly"] + "/outputs/scaffolds/"
    shell:
        """
        cat {input} | awk '{{if (substr($0, 1, 1)=='>') {{filename=(substr($0,2) '.fasta'}} print $0 >> {params}/filename
        close({params}/filename)}}'
        """"


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
# TODO
