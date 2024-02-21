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
        # report = config["assembly"] + "/reports/quast/report.html",
        report(
            config["assembly"] + "/reports/quast/report.html", 
            caption="../docs/captions/quast.rst", 
            category="Descriptive Stats")
    params:
        out_pfx = f"{config['assembly']}/reports/quast/",
        threads = config["threads"],
        log = f"{config['assembly']}/logs/quast.log",
        quast_path = config['quast_path']
    shell:
        """
        {params[quast_path]}/quast.py --large {input[assembly]} --glimmer -b --threads {params[threads]} -L --pacbio {input[reads]} -o {params[out_pfx]} 2> {params[log]} &&
        cp {params[out_pfx]}quast.log {params[log]}
        """


        
