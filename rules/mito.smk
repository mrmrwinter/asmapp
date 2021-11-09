##############################################################################

# MITOCHONDRIAL DETECTION

# rule mito_identification:
#     input:
#         initial = "data/assemblies/" + config["assembly"] + ".fasta",
#           initial_db =
#   output:
#       config["assembly"] + "reports/blast/mito_blast.out"
# params:
#     out_pfx = config["assembly"] + "reports/blast/
#     mito_ref = "data/"
#       threads =
# shell:
        # "blastn -query {params[mito_ref]} -db {input[initial_db]} -outfmt 6 -max_target_seqs 1 -out {params[out_pfx]} + "mito_blast.out -num_threads {params[threads]}"




rule mito_tagging:
    input:
        initial = "data/assemblies/" + config["assembly"] + ".fasta",
        config["assembly"] + "reports/blast/mito_blast.out"
    output:
        mito_tagged = "data/assemblies/" + config["assembly"] + ".mito_tagged.fasta",
        no_mito = "data/assemblies/" + config["assembly"] + ".no_mito.fasta",
    params:
        mito =
    run:
        # read in the header of the mito scaffold
        import pandas as import pd

        pd.read_csv(input[1], sep = "\t")

        mito_tig = value[1]

        # code to tag the mito contig header
        from Bio import SeqIO
        to_add = "mitochondrial_"
        with open(output[0], "w") as outputs:
            for r in SeqIO.parse(input[0], "fasta"):
                if r.id = str(mito_tig)
                r.id = (to_add + r.description).replace(" ", "_")
                r.description = r.id
                SeqIO.write(r, outputs, "fasta")

        # code to remove the mito contig from the assembly
        
