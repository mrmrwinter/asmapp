
##############################################################################

# CEGMA

rule CEGMA:
    container:
        "docker://chrishah/cegma"
    input:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",
    output:
        config["assembly"] + "/outputs/cegma/" + config["assembly"] + ".completeness_report"
    params:
        config["assembly"] + "/outputs/cegma/" + config["assembly"],
        threads = config["threads"]
    shell:
        """
        export CEGMATMP='{params[0]}'
        cegma --threads {params[threads]} --genome {input[0]} -o {params[0]}
        """
