##############################################################################

# PAIRS TABLE with BLAST

rule make_blast_database_initial:  # Rule to make database of cds fasta
    input:
        "data/assemblies/" + config["assembly"] + ".fasta" # input to the rule
    output:
        nhr = "data/databases/" + config["assembly"] + "/initial_" + config["assembly"] + ".nhr",   # all outputs expected from the rule
        nin = "data/databases/" + config["assembly"] + "/initial_" + config["assembly"] + ".nin",
        nsq = "data/databases/" + config["assembly"] + "/initial_" + config["assembly"] + ".nsq"
    params:
        "data/databases/" + config["assembly"] + "/initial_" + config["assembly"]   # prefix for the outputs, required by the command
    shell:  # shell command for the rule
        "makeblastdb \
        -in {input} \
        -out {params} \
        -dbtype nucl"  # the database type


rule blast_nonself_initial:
    input:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",
        db = "data/databases/" + config["assembly"] + "/initial_" + config["assembly"] + ".nin"
    output:
        tsv = config["assembly"] + "/outputs/blast/initial_blast.out"
    params:
        out_pfx = config["assembly"] + "/outputs/blast",
        db_pfx = "data/databases/" + config["assembly"] + "/initial_" + config["assembly"],
        threads = config["threads"]
    shell:
        "blastn -query {input[0]} -db {params[1]} -outfmt 6 -max_target_seqs 2 -out {params[0]}/initial_blast.out -num_threads {params[threads]}"


##############################################################################
