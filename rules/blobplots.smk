# BLOBPLOTS

# Blast the contigs against the NCBI nt database to assign taxonomy to them
rule tax_blast:
    input:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",
        # os.path.join(config["ncbi_nt_path"], "nt.115.nin")
    output:
        config["assembly"] + "/reports/blast/contaminant_taxonomy.blast.out"
    params:
        threads = config["threads"],
        log = f"{config['assembly']}/logs/tax_blast.log",
        ncbi_db = config['ncbi_db'],
        ncbi_db_path = config['ncbi_db_path']
    shell:
        """
        export BLASTDB={params[ncbi_db_path]}

        blastn \
        -db {params[ncbi_db]} \
        -query {input[assembly]} \
        -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -evalue 1e-25 \
        -num_threads {params[threads]} \
        -out {output} 2> {params[log]}
        """

# Create the intial blobject using the mapped reads, the contig taxonomy results, and the assembly.
rule blob_create:
    input:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",
        reads = f"{config['assembly']}/outputs/mapping/{config['reads']}.sorted.bam",
        hits = f"{config['assembly']}/reports/blast/contaminant_taxonomy.blast.out",
        index = f"{config['assembly']}/outputs/mapping/{config['reads']}.sorted.bam.bai"
    output:
        config["assembly"] + "/reports/blobtools/" + config["assembly"] + ".blobDB.json"
    params:
        out = config["assembly"] + "/reports/blobtools/" + config["assembly"],
        db = config["ncbi_db_path"],
        log = f"{config['assembly']}/logs/blob_create.log",
    shell:
        "blobtools create \
        -i {input.assembly} \
        -b {input.reads} \
        -t {input.hits} \
        -o {params[out]} \
        --names {params[db]}/names.dmp \
        --nodes {params[db]}/nodes.dmp 2> {params[log]}"


# View the blobject as a table
rule blobtools_view:
    input:
        blob_json = config["assembly"] + "/reports/blobtools/" + config["assembly"] + ".blobDB.json"
    output:
        report(
            config["assembly"] + "/reports/blobtools/" + config["assembly"] + ".blobDB.table.txt",
            caption="../docs/captions/blobtools.rst",
            category="Contamination reports"
        )
    params:
        out_pfx = config["assembly"] + "/reports/blobtools/",
        log = f"{config['assembly']}/logs/blobtools_view.log",
    shell:
        "blobtools view \
        -i {input[blob_json]} \
        --out {params[out_pfx]} 2> {params[log]}"


# Plot the blobplot
rule blobtools_plot:
    input:
        blob_json = config["assembly"] + "/reports/blobtools/" + config["assembly"] + ".blobDB.json"
    output:
        report(
            config["assembly"] + "/reports/blobtools/" + config["assembly"] + ".blobDB.json.bestsum.phylum.p8.span.100.blobplot.bam0.png", 
            caption="../docs/captions/blobtools.rst", 
            category="Contamination reports")
    params:
        out = config["assembly"] + "/reports/blobtools/",
        log = f"{config['assembly']}/logs/blobtools_plot.log",
    shell:
        "blobtools plot \
        -i {input[blob_json]} \
        --out {params[out]} 2> {params[log]}"
