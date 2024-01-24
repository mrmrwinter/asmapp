# COMPLETENESS ANALYSIS

# Run CEGMA on the assembly
rule CEGMA:
    container:
        "docker://chrishah/cegma"
    input:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",
    output:
        # config["assembly"] + "/reports/cegma/" + config["assembly"] + ".completeness_report"
        report(
            config["assembly"] + "/reports/cegma/" + config["assembly"] + ".completeness_report",
            caption="../reports/CEGMA.rst",
            category="Completeness"
        )
    params:
        out_pfx = config["assembly"] + "/reports/cegma/" + config["assembly"],
        threads = config["threads"]
    shell:
        """
        export CEGMATMP='{params[out_pfx]}'
        cegma --threads {params[threads]} --genome {input[assembly]} -o {params[out_pfx]}
        """
