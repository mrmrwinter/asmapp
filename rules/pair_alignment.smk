# PAIRS ALIGNMENTS FROM NUCMER

rule only_pairs:
    input:
       blast = config["assembly"] + "/reports/blast/blast.out",
    output:
       only_pairs_table = config["assembly"] + "/reports/blast/blast.onlyPairs.tsv"
    run:
       import pandas as pd

       blast_output = pd.read_csv(input[0], sep="\t", header = None) # snakemake.input[0] is the blast table

       initial_pairs = pd.DataFrame(columns = ['query', 'hit'])

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

#
rule single_scaffold_extraction:
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta",
    output:
        directory(config["assembly"] + "/tmp")
    params:
        out_dir = config["assembly"] + "/tmp/"
    run:
        import os

        os.system("mkdir -p " + params[0])

        os.system("awk '/^>/ {OUT='" + params[0] + "' substr($0,2) '.fasta'}; OUT{print >OUT}' " + input[0])


#
#
# # nucmer analysis/alignment
rule nucmer_alignment:
    input:
        directory(config["assembly"] + "tmp/"),
        only_pairs_table = config["assembly"] + "/reports/blast/blast.onlyPairs.tsv"
    output:
        out_dir = directory(config["assembly"] + "reports/nucmer/")
    params:
        tigs = config["assembly"] + "tmp/",
        out_dir = config["assembly"] + "reports/nucmer/"
    run:
        import glob
        import os

        # os.system("mkdir nucmer_purged/")

        pairs = pd.read_csv()

        for index, value in only_pairs.iterrows():
            q = params[0] + value[0] + ".fasta"
            h = params[0] + value[1] + ".fasta"
            nucmer = "nucmer -p " + params[1] + "nucmer_initial/nucmer." + str(index) + " " + q + " " + h
            os.system(nucmer)


        for delta in glob.glob("nucmer_initial/*.delta"):
            pair = delta.replace("nucmer_initial/","").replace(".delta","")
            mummer = "mummerplot -l -f --png --large " + delta + " -p " + params[0] + "nucmer_initial/" + pair
            os.system(mummer)



# initial
# PAIRS ALIGNMENTS FROM NUCMER

rule only_pairs_initial:
    input:
       blast = config["assembly"] + "/reports/blast/initial_blast.out",
    output:
       only_pairs_table = config["assembly"] + "/reports/blast/initial_blast.onlyPairs.tsv"
    run:
       import pandas as pd

       blast_output = pd.read_csv(input[0], sep="\t", header = None) # snakemake.input[0] is the blast table

       initial_pairs = pd.DataFrame(columns = ['query', 'hit'])

       for index, value in blast_output.iterrows():
           if value[0] != value[1]:
               initial_pairs.loc[index, ['query']] = value[0]
               initial_pairs.loc[index, ['hit']] = value[1]

       only_initial_pairs = initial_pairs.drop_duplicates()

       only_initial_pairs.to_csv(output[0], sep='\t')

# PRETTY SURE THIS LOOP IS INEFFICIENT AS ...


# IS THERE A WAY TO TAKE A SAMPLE LIST FROM THE ONLY PAIRS COLUMN HERE?
# IF SO I WONT HAVE TO USE THE DIRECTORY OUTPUT FLAG WHICH IS A pain
# CAN INSTEAD LOOP THROUGH THE LIST FOR EACH INDIVIDUAL ASSEMBLY AS THE LIST OF PAIRS WILL BE NOVEL AND INDIVIDUAL

#
rule single_scaffold_extraction:
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta",
    output:
        directory(config["assembly"] + "/tmp_initial")
    params:
        out_dir = config["assembly"] + "/tmp_initial/"
    run:
        import os

        os.system("mkdir -p " + params[0])

        os.system("awk '/^>/ {OUT='" + params[0] + "' substr($0,2) '.fasta'}; OUT{print >OUT}' " + input[0])


#
#
# # nucmer analysis/alignment
rule nucmer_alignment:
    input:
        directory(config["assembly"] + "tmp/"),
        only_pairs_table = config["assembly"] + "/reports/blast/blast.onlyPairs.tsv"
    output:
        out_dir = config["assembly"] + "reports/nucmer/"
    params:
        tigs = config["assembly"] + "tmp/",
        out_dir = config["assembly"] + "reports/nucmer/"
    run:
        import glob
        import os

        # os.system("mkdir nucmer_initial_purged/")

        pairs = pd.read_csv()

        for index, value in only_pairs.iterrows():
            q = params[0] + value[0] + ".fasta"
            h = params[0] + value[1] + ".fasta"
            nucmer = "nucmer -p " + params[1] + "nucmer_initial/nucmer." + str(index) + " " + q + " " + h
            os.system(nucmer)


        for delta in glob.glob("nucmer_initial/*.delta"):
            pair = delta.replace("nucmer_initial/","").replace(".delta","")
            mummer = "mummerplot -l -f --png --large " + delta + " -p " + params[0] + "nucmer_initial/" + pair
            os.system(mummer)
