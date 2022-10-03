##############################################################################

# MITOCHONDRIAL DETECTION

rule mito_identification:
    input:
        db = "data/databases/" + config["assembly"] + "/" + config["assembly"] + ".nin",
        mito_ref = "data/assemblies/mitoref.fasta"
    output:
        config["assembly"] + "reports/blast/mito_blast.out"
    params:
        db = "data/databases/" + config["assembly"] + "/" + config["assembly"],
        out_pfx = config["assembly"] + "reports/blast/",
        threads = config["threads"]
    shell:
        "blastn -query {input[0]} -db {params[db]} -outfmt 6 -max_target_seqs 1 -out {params[out_pfx]} =o {output} -num_threads {params[threads]}"


rule mito_tagging:
    input:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",
        mito_blast = config["assembly"] + "reports/blast/mito_blast.out"
    output:
        mito_tagged = config["assembly"] + "/outputs/assemblies/" + config["assembly"] + ".mito_tagged.fasta",
        no_mito = config["assembly"] + "/outputs/assemblies/" + config["assembly"] + ".no_mito.fasta",
    run:
        import pandas as pd
        from Bio import SeqIO

        # read in the header of the mito scaffold
        blast_output = pd.read_csv(input[1], sep="\t", header = None)
        for index, value in blast_output.iterrows():
            if index == 0:
                mito_tig = value[1]

        # code to tag the mito contig header
        to_add = "mitochondrial_"

        with open(output[0], "w") as mito_tagged:
            for r in SeqIO.parse(input[0], "fasta"):
                if r.id in str(mito_tig):
                    r.id = (to_add + r.description).replace(" ", "_")
                    r.description = r.id
                    SeqIO.write(r, mito_tagged, "fasta")

        # code to remove the mito contig from the assembly
        with open(output[1], "w") as no_mito: # etc
            for seq in SeqIO.parse(input[0], "fasta"):
                if seq.id not in str(mito_tig):
                # write to new file
                    SeqIO.write(seq, no_mito, "fasta")
        