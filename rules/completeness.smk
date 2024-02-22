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
        out_pfx = config["assembly"] + "/reports/cegma/" + config["assembly"],
        threads = config["threads"],
        log = f"{config['assembly']}/logs/CEGMA.log",
    shell:
        """
        export CEGMATMP='{params[out_pfx]}'
        cegma --threads {params[threads]} --genome {input[assembly]} -o {params[out_pfx]} 2> {params[log]}
        """


# Run CEGMA on the assembly
# Run CEGMA on the assembly
rule BUSCO:
    conda:
        "../envs/BUSCO.yaml"
    input:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",
    output:
        # config["assembly"] + "/reports/cegma/" + config["assembly"] + ".completeness_report"
        report(
            f"{config['assembly']}/reports/busco/run_{config['busco_lineage']}/short_summary.txt",
            caption="../docs/captions/BUSCO.rst",
            category="Completeness"
        )
    params:
        out_pfx = config["assembly"] + "/reports/busco/",
        threads = config["threads"],
        log = f"{config['assembly']}/logs/BUSCO.log",
        busco_lineage = config['busco_lineage']
    shell:
        """
        busco -i {input[assembly]} -l {params[busco_lineage]} -o {params[out_pfx]} -m genome --cpu {params[threads]} -f 2> {params[log]}
        mv busco*.log {params[out_pfx]}
        """
