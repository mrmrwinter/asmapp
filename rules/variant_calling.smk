# VARIANT CALLING

# Run sniffles on the assembly to detect structural variation
rule sniffles:
    input:
        initial="data/assemblies/" + config["assembly"] + ".fasta",
        bam=config["assembly"] + "/outputs/initial/initial_asm.sorted.bam",
    output:
        config["assembly"] + "/outputs/variant_calling/" + config["assembly"] + ".vcf",
    params:
        threads=config["threads"],
    shell:
        "sniffles -i {input[1]}  --vcf {output} --reference {input[0]}"