# NUCMER ASSEMBLY DOTPLOTS


# Perform a self-by-self alignment with nucmer and print a dotplot
rule nucmer_self:
    input:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta"
    output:
        report(
            config["assembly"] + "/outputs/nucmer/nucmer.self_v_self.png",
            caption = "../docs/captions/nucmer.rst", 
            category = "Dotplots"
            ),
        config["assembly"] + "/outputs/nucmer/nucmer.self_v_self.delta"
    params:
        out_pfx = "nucmer.self_v_self",
        out_path = config["assembly"] + "/outputs/nucmer/",
        log = f"{config['assembly']}/logs/{rule}.log",
    shell:
        """
        mkdir -p tmp/
        nucmer -p tmp/{params[out_pfx]} {input[assembly]} {input[assembly]}
        cp tmp/{params[out_pfx]}.delta {params[out_path]}
        mummerplot -l -f --png --large {params[out_path]}{params[out_pfx]}.delta -p {params[out_path]}{params[out_pfx]} 2> {params[log]}
        """


# Perform a self-by-reference alignment with nucmer and print a dotplot
rule nucmer_initial_vs_reference:
    input:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",
        reference = "data/assemblies/" + config["reference"] + ".fasta",
        # nuc_fai = config["assembly"] + "/outputs/nucmer/scaffolds.reduced.fasta.fai"
    output:
        report(
            config["assembly"] + "/outputs/nucmer/nucmer.self_v_ref.png", 
            caption = "../docs/captions/nucmer.rst", 
            category = "Dotplots"
            ),
        config["assembly"] + "/outputs/nucmer/nucmer.self_v_ref.delta"
    params:
        out_pfx = "nucmer.self_v_ref",
        out_path = config["assembly"] + "/outputs/nucmer/",
        reference = config["reference"] + ".fasta",
        log = f"{config['assembly']}/logs/{rule}.log",
    shell:
        """
        mkdir -p tmp/
        nucmer -p tmp/{params[out_pfx]} {input[assembly]} {input[reference]}
        cp tmp/{params[out_pfx]}.delta {params[out_path]}
        mummerplot -l -f --png --large {params[out_path]}{params[out_pfx]}.delta -p {params[out_path]}{params[out_pfx]} 2> {params[log]}
        """

