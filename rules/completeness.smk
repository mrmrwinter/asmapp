# COMPLETENESS ANALYSIS


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


# rule busco:
#     conda:
#         "../envs/BUSCO.yaml"
#     input:
#         assembly = "data/assemblies/" + config["assembly"] + ".fasta",
#     output:
#         report(
#             config["assembly"] + "/reports/BUSCO/run_{params[lineage]}/short_summary." + config["assembly"] + ".txt",
#             caption="../reports/BUSCO.rst",
#             category="Completeness"
#         )
#     params:
#         lineage = config["busco_lineage"],
#         out_pfx = config["assembly"] + "/reports/BUSCO/"
#     shell:
#         "busco -m genome -i {input} -o {params[out_pfx]} -l {params}"
