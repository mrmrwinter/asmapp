# NUCMER ASSEMBLY DOTPLOTS

# Perform a self-by-self alignment with nucmer and print a dotplot
rule nucmer_self:
    input:
        initial = "data/assemblies/" + config["assembly"] + ".fasta"
    output:
        report(
            config["assembly"] + "/outputs/nucmer/nucmer.self_v_self.png",
            caption = "reports/nucmer.rst", 
            category = "Dotplots"
            ),
        config["assembly"] + "/outputs/nucmer/nucmer.self_v_self.delta"
    params:
        "nucmer.self_v_self",
        config["assembly"] + "/outputs/nucmer/",
    shell:
        """
        mkdir -p tmp/
        nucmer -p tmp/{params[0]} {input[0]} {input[0]}
        cp tmp/{params[0]}.delta {params[1]}
        mummerplot -l -f --png --large {params[1]}{params[0]}.delta -p {params[1]}{params[0]}
        """

# Perform a self-by-reference alignment with nucmer and print a dotplot
rule nucmer_initial_vs_reference:
    input:
        initial = "data/assemblies/" + config["assembly"] + ".fasta",
        reference = config["reference"] + ".fasta.gz",
        # nuc_fai = config["assembly"] + "/outputs/nucmer/scaffolds.reduced.fasta.fai"
    output:
        report(
            config["assembly"] + "/outputs/nucmer/nucmer.self_v_ref.png", 
            caption = "reports/nucmer.rst", 
            category = "Dotplots"
            ),
        config["assembly"] + "/outputs/nucmer/nucmer.self_v_ref.delta"
    params:
        "nucmer.self_v_ref",
        config["assembly"] + "/outputs/nucmer/",
        config["reference"] + ".fasta"
    shell:
        """
        mkdir -p tmp/
        gunzip -c {input[reference]} > {params[2]}
        nucmer -p tmp/{params[0]} {input[0]} {params[2]}
        cp tmp/{params[0]}.delta {params[1]}
        mummerplot -l -f --png --large {params[1]}{params[0]}.delta -p {params[1]}{params[0]}
        """

