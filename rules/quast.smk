# QUAST

 # Performing QUAST assembly appraisal
rule quast:
    message:
        "[INFO] Performing QUAST appraisal on assemblies..."
    input:
        initial = "data/assemblies/" + config["assembly"] + ".fasta",
        # assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fa",
        reads = "data/reads/" + config["reads"] + ".fasta",
        reference = config["reference"] + ".fasta.gz"
    output:
        # config["assembly"] + "/reports/quast/report.html",
        report(
            config["assembly"] + "/reports/quast/report.html", 
            caption="../reports/quast.rst", 
            category="Descriptive Stats")
    params:
        out_pfx = config["assembly"] + "/reports/quast/",
        threads = config["threads"]
    shell:
        config["quast_path"] + "/quast.py --large {input[0]} --glimmer -b --threads {params[1]} -L --pacbio {input[reads]} -o {params[0]}"


        
