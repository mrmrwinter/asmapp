# QUAST

# Performing QUAST assembly appraisal
rule quast:
    message:
        "[INFO] Performing QUAST appraisal on assemblies..."
    input:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",
        # assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fa",
        reads = "data/reads/" + config["reads"] + ".fasta",
        reference = "data/assemblies/" + config["reference"] + ".fasta"
    output:
        # config["assembly"] + "/reports/quast/report.html",
        report(
            config["assembly"] + "/reports/quast/report.html", 
            caption="../docs/captions/quast.rst", 
            category="Descriptive Stats")
    params:
        out_pfx = config["assembly"] + "/reports/quast/",
        threads = config["threads"],
        log = f"{config['assembly']}/logs/{rule}.log",
    shell:
        "quast --large {input[assembly]} --glimmer -b --threads {params[threads]} -L --pacbio {input[reads]} -o {params[out_pfx]} 2> {params[log]}"


        
