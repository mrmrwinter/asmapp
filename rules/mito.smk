# MITOCHONDRIAL DETECTION


# Identify if the mitochondrial genome is present in the assembly
rule mito_identification:
    input:
        nin = f"data/databases/{config['assembly']}/{config['assembly']}.nin",
        mito_ref = f"data/assemblies/{config['mitochondrial']}.fasta"
    output:
        f"{config['assembly']}/reports/blast/mito_blast.out"
    params:
        db = "data/databases/" + config["assembly"] + "/" + config["assembly"],
        out_pfx = config["assembly"] + "/reports/blast/",
        threads = config["threads"],
        log = f"{config['assembly']}/logs/mito_identification.log",
    shell:
        "blastn -query {input[mito_ref]} -db {params[db]} -outfmt 6 -max_target_seqs 1 -out {output} -num_threads {params[threads]} 2> {params[log]}"


# If found, tag the mitochondrial genome
rule mito_tagging:
    input:
        assembly = f"data/assemblies/{config['assembly']}.fasta",
        mito_blast = f"{config['assembly']}/reports/blast/mito_blast.out"
    output:
        mito_tagged = config["assembly"] + "/outputs/assemblies/" + config["assembly"] + ".mito_tagged.fasta",
        no_mito = config["assembly"] + "/outputs/assemblies/" + config["assembly"] + ".no_mito.fasta",
    run:
        import traceback
        
        try:
            import pandas as pd
        except ImportError as e:
            print(f"Error importing pandas: {e}")
            traceback.print_exc()
            
        try:
            from Bio import SeqIO
        except ImportError as e:
            print(f"Error importing Bio.SeqIO: {e}")
            traceback.print_exc()

        # read in the header of the mito scaffold
        blast_output = pd.read_csv(f"{input['mito_blast']}", sep="\t", header = None)

        for index, value in blast_output.iterrows():
            if index != 0:
                continue
            else:
                mito_tig = value[1]
                break
        # code to tag the mito contig header
        to_add = "mitochondrial_"

        with open(f"{output['mito_tagged']}", "w") as mito_tagged:
            for r in SeqIO.parse(input['assembly'], "fasta"):
                print(r.id)
                if r.id == str(mito_tig):
                    r.id = (f"{to_add}{r.id}").replace(" ", "_")
                    SeqIO.write(r, mito_tagged, "fasta")
        print("Mitochondrial scaffold succesfully tagged.")
        # code to remove the mito contig from the assembly
        with open(f"{output['no_mito']}", "w") as no_mito: # etc
            for seq in SeqIO.parse(input[0], "fasta"):
                if seq.id not in str(mito_tig):
                # write to new file
                    SeqIO.write(seq, no_mito, "fasta")
        print("Mito-free assembly generated.")
       