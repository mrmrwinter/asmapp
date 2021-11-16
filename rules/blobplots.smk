
##############################################################################

# BLOBPLOTS


# setup
# requires that the ncbi database is installed and configured in the config files
# needs taxdumpdownloading and unpacking in the same folder as the nt database

rule blast:
    input:
        "data/assemblies/" + config["assembly"] + ".fasta",
    output:
        config["assembly"] + "/reports/blast/contaminant_taxonomy.blast.out"
    params:
        blastdb=config["ncbi_nt_path"],
        threads = config["threads"]
    shell:
        "blastn \
        -db {params.blastdb}/nt \
        -query {input} \
        -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -evalue 1e-25 \
        -num_threads {params[threads]} \
        -out {output}"


rule blob_create:
    conda:
        "../envs/blobtools.yaml"
    input:
        initial = "data/assemblies/" + config["assembly"] + ".fasta",
        reads = config["assembly"] + "/outputs/initial/initial_asm.sorted.bam",
        hits = config["assembly"] + "/reports/blast/contaminant_taxonomy.blast.out"
    output:
        config["assembly"] + "/reports/blobtools/" + config["assembly"] + ".blobDB.json"
    params:
        out = config["assembly"] + "/reports/blobtools/" + config["assembly"],
        db = config["ncbi_nt_path"],
    shell:
        "blobtools create \
        -i {input.initial} \
        -b {input.reads} \
        -t {input.hits} \
        -o {params.out} \
        --names {params[1]}/names.dmp \
        --nodes {params[1]}/nodes.dmp"


rule blobtools_view:
    conda:
        "../envs/blobtools.yaml"
    input:
        config["assembly"] + "/reports/blobtools/" + config["assembly"] + ".blobDB.json"
    output:
        config["assembly"] + "/reports/blobtools/" + config["assembly"] + ".blobDB.table.txt"
    params:
        config["assembly"] + "/reports/blobtools/" + config["assembly"]
    shell:
        "blobtools view \
        -i {input} \
        --out {params[0]}"


rule blobtools_plot:
    conda:
        "../envs/blobtools.yaml"
    input:
        config["assembly"] + "/reports/blobtools/" + config["assembly"] + ".blobDB.json"
    output:
        config["assembly"] + "/reports/blobtools/" + config["assembly"] + ".blobDB.json.bestsum.phylum.p8.span.100.blobplot.bam0.png"
    params:
        out = config["assembly"] + "/reports/blobtools/"
    shell:
        "blobtools plot \
        -i {input} \
        --out {params[out]}"
