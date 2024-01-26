# VARIANT CALLING

# Run sniffles on the assembly to detect structural variation
rule sniffles:
    input:
        bai = f"{config['assembly']}/outputs/mapping/{config['reads']}.sorted.bam.bai",
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",
        bam = f"{config['assembly']}/outputs/mapping/{config['reads']}.sorted.bam",
    output:
        vcf = f"{config['assembly']}/outputs/variant_calling/{config['assembly']}_{config['reads']}.vcf",
    params:
        threads = config["threads"],
    shell:
        "sniffles -i {input[bam]}  --vcf {output} --reference {input[assembly]}"