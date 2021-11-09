
##############################################################################

# BLOBPLOTS


#####   BLASTing   #####

configfile: "config.yaml"

# rule all:
#     input:
        config["assembly"] + "reports/blast/contaminant_taxonomy.blast.out"
        config["assembly"] + "reports/blobtools/" + config["assembly"] + "blobDB.json"




rule blast:
    conda:
        "envs/blast.yaml"
    input:
        "data/assemblies/" + config["assembly"] + ".fasta",
    output:
        config["assembly"] + "reports/blast/contaminant_taxonomy.blast.out"
    threads:
        congi["threads"]
    params:
        blastdb="nt"
    shell:
        "blastn \
        -db {params.blastdb} \
        -query {input} \
        -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -evalue 1e-25 \
        -num_threads {threads} \
        -out {output}"


rule blob_create:
    conda:
        "envs/blobtools.yaml"
    input:
        initial = "data/assemblies/" + config["assembly"] + ".fasta",
        reads = "data/reads/" + config["reads"] + ".fasta",
        hits = config["assembly"] + "reports/blast/contaminant_taxonomy.blast.out"
    output:
        config["assembly"] + "reports/blobtools/" + config["assembly"] + "blobDB.json"
    params:
        out = config["assembly"] + "reports/blobtools/" + config["assembly"]
    threads:
        28
    shell:
        config["blobtools_path"] + "blobtools create \
        -i {input.initial} \
        -b {input.reads} \
        -t {input.hits} \
        -o {params.out}"


rule blobtools_view:
    conda:
        "envs/blobtools.yaml"
    input:
        config["assembly"] + "reports/blobtools/" + config["assembly"] + ".blobDB.json"
    output:
        config["assembly"] + "reports/blobtools/" + config["assembly"] + ".blobDB.table.txt"
    params:
        config["assembly"] + "reports/blobtools/" + config["assembly"]
    threads:
        28
    shell:
        config["blobtools_path"] + "blobtools view \
        -i {input} \
        --out {params[0]}"


rule blobtools_plot:
    conda:
        "envs/blobtools.yaml"
    input:
        config["assembly"] + "reports/blobtools/" + config["assembly"] + ".blobDB.json"
    output:
        config["assembly"] + "reports/blobtools/" + config["assembly"] + ".blobDB.json.bestsum.phylum.p8.span.100.blobplot.bam0.png"
    params:
        out="reports/blobtools/" + config["assembly"] + "."
    threads:
        28
    shell:
        config["blobtools_path"] + "blobtools plot \
        -i {input} \
        --out {params.out}"
