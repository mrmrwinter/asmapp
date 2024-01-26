# BLOBPLOTS

# Blast the contigs against the NCBI nt database to assign taxonomy to them
rule tax_blast:
    input:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",
        # os.path.join(config["ncbi_nt_path"], "nt.115.nin")
    output:
        config["assembly"] + "/reports/blast/contaminant_taxonomy.blast.out"
    params:
        blastdb = config["ncbi_nt_path"],
        threads = config["threads"]
    shell:
        "blastn \
        -db {params[blastdb]} \
        -query {input[assembly]} \
        -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -evalue 1e-25 \
        -num_threads {params[threads]} \
        -out {output}"

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
        db = config["ncbi_nt_path"],
    shell:
        "blobtools create \
        -i {input.assembly} \
        -b {input.reads} \
        -t {input.hits} \
        -o {params[out]} \
        --names {params[db]}/names.dmp \
        --nodes {params[db]}/nodes.dmp"


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
        out_pfx = config["assembly"] + "/reports/blobtools/"
    shell:
        "blobtools view \
        -i {input[blob_json]} \
        --out {params[out_pfx]}"


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
        out = config["assembly"] + "/reports/blobtools/"
    shell:
        "blobtools plot \
        -i {input[blob_json]} \
        --out {params[out]}"

