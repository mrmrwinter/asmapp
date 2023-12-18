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
            caption="../docs/captions/CEGMA.rst",
            category="Completeness"
        )
    params:
        config["assembly"] + "/reports/cegma/" + config["assembly"],
        threads = config["threads"]
    shell:
        """
        export CEGMATMP='{params[0]}'
        cegma --threads {params[threads]} --genome {input[0]} -o {params[0]}
        """
