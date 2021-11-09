##############################################################################

# MITOCHONDRIAL DETECTION

# rule mito_identification:
#     input:
#         initial = "data/assemblies/" + config["assembly"] + ".fasta",
#
#



rule mito_tagging:
    input:

    output:

    run:
        from Bio import SeqIO
        to_add = "mitochondrial_"
        with open(output[0], "w") as outputs:
            for r in SeqIO.parse(input[0], "fasta"):
                if r.id = "xxxxxxx"
                r.id = (to_add + r.description).replace(" ", "_")
                r.description = r.id
                SeqIO.write(r, outputs, "fasta")
