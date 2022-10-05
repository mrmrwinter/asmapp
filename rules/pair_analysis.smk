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
        tsv = config["assembly"] + "/reports/blast/blast.out"
    params:
        out_pfx = config["assembly"] + "/reports/blast",
        db_pfx = "data/databases/" + config["assembly"] + "/" + config["assembly"],
        threads = config["threads"]
    shell:
        "blastn -query {input[0]} -db {params[1]} -outfmt 6 -max_target_seqs 2 -out {params[0]}/blast.out -num_threads {params[threads]}"



rule only_pairs:
    input:
       blast = config["assembly"] + "/reports/blast/blast.out",
    output:
       only_pairs_table = config["assembly"] + "/reports/pairs_analysis/blast/blast.onlyPairs.tsv"
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

# PAIRS WITH MASH

rule mash_sketch:
    input:
        scaffolds = config["assembly"] + "/outputs/scaffolds/{all_scaffs}.fasta",
    output:
        config["assembly"] + "reports/pair-analysis/mash/" + config["assembly"] + ".sketch",
    params:
        scaffolds = config["assembly"] + "/outputs/scaffolds/",
    shell:
        "mash sketch {params[scaffolds]}/*.fasta -o {output}"


rule mash_dist:
    input:
        config["assembly"] + "reports/pair-analysis/mash/" + config["assembly"] + ".sketch",
    output:
        config["assembly"] + "reports/pair-analysis/mash/" + config["assembly"] + ".dist",
    params:
        scaffolds = config["assembly"] + "/outputs/scaffolds/",
    shell:
        "mash sketch {params}/*.fasta -o {output}"

rule get_mash_distances:
    input:
    output:
    shell:
    









##################################

# # nucmer analysis/alignment
rule nucmer_alignment:
    input:
        # directory(config["assembly"] + "tmp_initial/"),
        only_pairs_table = config["assembly"] + "/reports/pairs_analysis/blast/blast.onlyPairs.tsv"
    output:
        report(
            directory(config["assembly"] + "reports/pair-analysis/nucmer/pairs"),
            caption="../docs/captions/pair_dotplots.rst",
            category="Pair analysis"
        )
    params:
        tigs = config["assembly"] + "tmp/",
        out_dir = config["assembly"] + "reports/nucmer/pairs/"
    run:
        import glob
        import os

        only_pairs = pd.read_csv(snakemake.input[1])

        for index, value in only_pairs.iterrows():
            q = params[0] + value[0] + ".fasta"
            h = params[0] + value[1] + ".fasta"
            nucmer = "nucmer -p " + params[1] + "nucmer/nucmer." + str(q.replace('.fasta','') + h.replace('.fasta','')) + " " + q + " " + h
            os.system(nucmer)


        for delta in glob.glob("nucmer/*.delta"):
            pair = delta.replace("nucmer/","").replace(".delta","")
            mummer = "mummerplot -l -f --png --large " + delta + " -p " + params[0] + "nucmer/" + pair
            os.system(mummer)



#######################################################################################

# dnadiff

# os.system("mkdir dnadiff_initial_purged/")
#
# for index, value in only_initial_pairs.iterrows():
#     dnadiff = "dnadiff -p dnadiff_initial_purged/nucmer." + str(index) + " -d nucmer_initial_purged/nucmer." + str(index) + ".delta"
#     os.system(dnadiff)
#
#