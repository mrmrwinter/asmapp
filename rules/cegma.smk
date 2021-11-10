
##############################################################################

# CEGMA

rule CEGMA_collapsed:
    container:
        "chrishah/cegma:2.5"
    input:
        initial = "data/assemblies/" + config["assembly"] + ".fasta",
    output:
        config["assembly"] + "outputs/cegma/" + config["assembly"] + ".completeness_report"
    params:
        config["assembly"] + "outputs/cegma/" + config["assembly"]
    shell:
        "cegma --genome {input[0]} -o {params[0]}"
