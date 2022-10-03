# NUCMER ASSEMBLY DOTPLOTS

rule nucmer_self:
    input:
        initial = "data/assemblies/" + config["assembly"] + ".fasta"
    output:
        report(config["assembly"] + "/outputs/nucmer/nucmer.self_v_self.png", caption = "outputs/nucmer_initial.rst", category = "Dotplots"),
        config["assembly"] + "/outputs/nucmer/nucmer.self_v_self.delta"
    params:
        "nucmer.initial",
        config["assembly"] + "/outputs/nucmer/",
    shell:
        """
        mkdir -p tmp/
        nucmer -p tmp/{params[0]} {input[0]} {input[0]}
        cp tmp/{params[0]}.delta {params[1]}
        mummerplot -l -f --png --large {params[1]}{params[0]}.delta -p {params[1]}{params[0]}
        """




rule nucmer_initial_vs_reference:
    input:
        initial = "data/assemblies/" + config["assembly"] + ".fasta",
        reference = config["reference"] + ".fasta.gz",
        # nuc_fai = config["assembly"] + "/outputs/nucmer/scaffolds.reduced.fasta.fai"
    output:
        config["assembly"] + "/outputs/nucmer/nucmer.self_v_ref.png",
        config["assembly"] + "/outputs/nucmer/nucmer.self_v_ref.delta"
    params:
        "nucmer.initial_v_ref",
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


 # TODO add nucmer delta transformation and synteny plotting

# rule nucmer_circles:
#     input:
#         fai = config["assembly"] + "/outputs/nucmer/scaffolds.reduced.fasta.fai"
#     output:
#         config["assembly"] + "/outputs/nucmer/circle.int_ref.png",
#     script:
#         "../scripts/snmk_nucmer_circles.R"

