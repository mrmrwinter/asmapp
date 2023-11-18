# VARIANT CALLING

# Run sniffles on the assembly to detect structural variation
rule sniffles:
    input:
        initial="data/assemblies/" + config['assembly'] + ".fasta",
        bam=f"{config['assembly']}/outputs/mapping/{config['reads']}.sorted.bam",
        bai = f"{config['assembly']}/outputs/mapping/{config['reads']}.sorted.bam.bai",
    output:
        vcf = f"{config['assembly']}/outputs/variant_calling/{config['assembly']}_{config['reads']}.vcf",
    params:
        threads=config["threads"],
    shell:
        "sniffles -i {input[1]}  --vcf {output} --reference {input[0]}"