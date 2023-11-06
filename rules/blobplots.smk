# BLOBPLOTS

# Blast the contigs against the NCBI nt database to assign taxonomy to them
rule tax_blast:
    input:
        "data/assemblies/" + config["assembly"] + ".fasta",
        # os.path.join(config["ncbi_nt_path"], "nt.115.nin")
    output:
        config["assembly"] + "/reports/blast/contaminant_taxonomy.blast.out"
    params:
        blastdb=config["ncbi_nt_path"],
        threads = config["threads"]
    shell:
        "blastn \
        -db {params.blastdb}nt \
        -query {input[0]} \
        -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -evalue 1e-25 \
        -num_threads {params[threads]} \
        -out {output}"

# Create the intial blobject using the mapped reads, the contig taxonomy results, and the assembly.
rule blob_create:
    input:
        initial = "data/assemblies/" + config["assembly"] + ".fasta",
        reads = f"{config['assembly']}/outputs/mapping/{config['reads']}.sorted.bam",
        hits = f"{config['assembly']}/reports/blast/contaminant_taxonomy.blast.out",
        index = f"{config['assembly']}/outputs/mapping/{config['reads']}.sorted.bam.bai"
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

# View the blobject as a table
rule blobtools_view:
    input:
        config["assembly"] + "/reports/blobtools/" + config["assembly"] + ".blobDB.json"
    output:
        report(
            config["assembly"] + "/reports/blobtools/" + config["assembly"] + ".blobDB.table.txt",
            caption="../docs/captions/blobtools.rst",
            category="Contamination reports"
        )
    params:
        config["assembly"] + "/reports/blobtools/"
    shell:
        "blobtools view \
        -i {input} \
        --out {params[0]}"

# Plot the blobplot
rule blobtools_plot:
    input:
        config["assembly"] + "/reports/blobtools/" + config["assembly"] + ".blobDB.json"
    output:
        report(
            config["assembly"] + "/reports/blobtools/" + config["assembly"] + ".blobDB.json.bestsum.phylum.p8.span.100.blobplot.bam0.png", 
            caption="../docs/captions/blobtools.rst", 
            category="Contamination reports")
    params:
        out = config["assembly"] + "/reports/blobtools/"
    shell:
        "blobtools plot \
        -i {input} \
        --out {params[out]}"