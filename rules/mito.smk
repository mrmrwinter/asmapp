# MITOCHONDRIAL DETECTION


# Identify if the mitochondrial genome is present in the assembly
rule mito_identification:
    input:
        nin = "data/databases/" + config["assembly"] + "/" + config["assembly"] + ".nin",
        mito_ref = config["mitochondrial"] + ".fasta"
    output:
        config["assembly"] + "/reports/blast/mito_blast.out"
    params:
        db = "data/databases/" + config["assembly"] + "/" + config["assembly"],
        out_pfx = config["assembly"] + "/reports/blast/",
        threads = config["threads"]
    shell:
        "blastn -query {input[mito_ref]} -db {params[db]} -outfmt 6 -max_target_seqs 1 -out {output} -num_threads {params[threads]}"


# If found, tag the mitochondrial genome
rule mito_tagging:
    input:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",
        mito_blast = config["assembly"] + "/reports/blast/mito_blast.out"
    output:
        mito_tagged = config["assembly"] + "/outputs/assemblies/" + config["assembly"] + ".mito_tagged.fasta",
        no_mito = config["assembly"] + "/outputs/assemblies/" + config["assembly"] + ".no_mito.fasta",
    run:
        import pandas as pd

        from Bio import SeqIO

        # read in the header of the mito scaffold
        blast_output = pd.read_csv(input['mito_blast'], sep="\t", header = None)

        for index, value in blast_output.iterrows():
            if index == 0:
                mito_tig = value[1]

        # code to tag the mito contig header
        to_add = "mitochondrial_"

        with open(output['mito_tagged'], "w") as mito_tagged:
            print("Mito_tagging output file created.")
            for r in SeqIO.parse(input['assembly'], "fasta"):
                if r.id == str(mito_tig):
                    print("Mito_tagging output file writing...")
                    r.id = (to_add + r.id).replace(" ", "_")
                    SeqIO.write(r, mito_tagged, "fasta")

        # code to remove the mito contig from the assembly
        with open(output['no_mito'], "w") as no_mito: # etc
            print("Mito purged assembly file created.")
            for seq in SeqIO.parse(input['assembly'], "fasta"):
                if seq.id not in str(mito_tig):
                    print("Mito purged assembly file writing...")
                # write to new file
                    SeqIO.write(seq, no_mito, "fasta")
