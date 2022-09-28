##############################################################################

# PAIRS TABLE with BLAST

rule make_blast_database:  # Rule to make database of cds fasta
    input:
        "data/assemblies/" + config["assembly"] + ".fasta" # input to the rule
    output:
        nhr = "data/databases/" + config["assembly"] + "/" + config["assembly"] + ".nhr",   # all outputs expected from the rule
        nin = "data/databases/" + config["assembly"] + "/" + config["assembly"] + ".nin",
        nsq = "data/databases/" + config["assembly"] + "/" + config["assembly"] + ".nsq"
    params:
        "data/databases/" + config["assembly"] + "/" + config["assembly"]   # prefix for the outputs, required by the command
    shell:  # shell command for the rule
        "makeblastdb \
        -in {input} \
        -out {params} \
        -dbtype nucl"  # the database type


rule blast_nonself:
    input:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",
        db = "data/databases/" + config["assembly"] + "/" + config["assembly"] + ".nin"
    output:
        tsv = config["assembly"] + "/outputs/blast/blast.out"
    params:
        out_pfx = config["assembly"] + "/outputs/blast",
        db_pfx = "data/databases/" + config["assembly"] + "/" + config["assembly"],
        threads = config["threads"]
    shell:
        "blastn -query {input[0]} -db {params[1]} -outfmt 6 -max_target_seqs 2 -out {params[0]}/blast.out -num_threads {params[threads]}"



rule only_pairs:
    input:
       blast = config["assembly"] + "/outputs/blast/blast.out",
    output:
       only_pairs_table = config["assembly"] + "/outputs/blast/blast.onlyPairs.tsv"
    run:
       import pandas as pd

       blast_output = pd.read_csv(input[0], sep="\t", header = None) # snakemake.input[0] is the blast table

       pairs = pd.DataFrame(columns = ['query', 'hit'])

       for index, value in blast_output.iterrows():
           if value[0] != value[1]:
               pairs.loc[index, ['query']] = value[0]
               pairs.loc[index, ['hit']] = value[1]

       only_pairs = pairs.drop_duplicates()

       only_pairs.to_csv(output[0], sep='\t')

# PRETTY SURE THIS LOOP IS INEFFICIENT AS ...

# IS THERE A WAY TO TAKE A SAMPLE LIST FROM THE ONLY PAIRS COLUMN HERE?
# IF SO I WONT HAVE TO USE THE DIRECTORY OUTPUT FLAG WHICH IS A pain
# CAN INSTEAD LOOP THROUGH THE LIST FOR EACH INDIVIDUAL ASSEMBLY AS THE LIST OF PAIRS WILL BE NOVEL AND INDIVIDUAL

##############################################################################
