# PAIRS TABLE with BLAST

# Rule to make database of assembly
rule make_blast_database:  
    input:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta" # input to the rule
    output:
        nhr = "data/databases/" + config["assembly"] + "/" + config["assembly"] + ".nhr",   # all outputs expected from the rule
        nin = "data/databases/" + config["assembly"] + "/" + config["assembly"] + ".nin",
        nsq = "data/databases/" + config["assembly"] + "/" + config["assembly"] + ".nsq"
    params:
        out_pfx = "data/databases/" + config["assembly"] + "/" + config["assembly"],
        log = f"{config['assembly']}/logs/make_blast_database.log",
    shell:  
        "makeblastdb \
        -in {input[assembly]} \
        -out {params[out_pfx]} \
        -dbtype nucl 2> {params[log]}"  # the database type


# BLAST the scaffolds back against the assembly
rule blast_nonself:
    input:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",
        db = "data/databases/" + config["assembly"] + "/" + config["assembly"] + ".nin"
    output:
        tsv = config["assembly"] + "/reports/blast/blast.out"
    params:
        out_pfx = config["assembly"] + "/reports/blast",
        db_pfx = "data/databases/" + config["assembly"] + "/" + config["assembly"],
        threads = config["threads"],
        log = f"{config['assembly']}/logs/blast_nonself.log",
    shell:
        "blastn -query {input[assembly]} -db {params[db_pfx]} -outfmt 6 -max_target_seqs 2 -out {params[out_pfx]}/blast.out -num_threads {params[threads]}"


# Remove self-to-self hits from the output
rule only_pairs:
    input:
       blast = config["assembly"] + "/reports/blast/blast.out",
    output:
       only_pairs_table = config["assembly"] + "/reports/pairs_analysis/blast/blast.onlyPairs.tsv"
    run:
       import pandas as pd

       blast_output = pd.read_csv(input[blast], sep="\t", header = None) 

       pairs = pd.DataFrame(columns = ['query', 'hit'])

       for index, value in blast_output.iterrows():
           if value[0] != value[1]:
               pairs.loc[index, ['query']] = value[0]
               pairs.loc[index, ['hit']] = value[1]

       only_pairs = pairs.drop_duplicates()

       only_pairs.to_csv(output[only_pairs_table], sep='\t')


# Perform nucmer alignment of potential pairs and print a dotplot for each
rule nucmer_pair_alignment:
    input:
        only_pairs_table = config["assembly"] + "/reports/pairs_analysis/blast/blast.onlyPairs.tsv"
    output:
        directory(config["assembly"] + "/reports/nucmer/pairs"),
        report(
            directory(config["assembly"] + "/reports/pairs_analysis/nucmer/pairs"),
            caption="../docs/captions/pair_dotplots.rst",
            category="Pair analysis"
        )
    params:
        tigs = config["assembly"] + "tmp/",
        out_dir = config["assembly"] + "/reports/nucmer/pairs/",
        log = f"{config['assembly']}/logs/nucmer_pair_alignment.log",
    run:
        import glob
        import os

        only_pairs = pd.read_csv(snakemake.input[only_pairs_table])

        for index, value in only_pairs.iterrows():
            q = params[tigs] + value[0] + ".fasta"
            h = params[tigs] + value[1] + ".fasta"
            nucmer = f"nucmer -p {params[out_dir]}nucmer/nucmer.{str(q.replace('.fasta','')}{h.replace('.fasta',''))} {q} {h} 2> {params[log]}"
            os.system(nucmer)


        for delta in glob.glob("nucmer/*.delta"):
            pair = delta.replace("nucmer/","").replace(".delta","")
            mummer = f"mummerplot -l -f --png --large {delta} -p {params[tigs]}nucmer/{pair} 2>> {params[log]}"
            os.system(mummer)

