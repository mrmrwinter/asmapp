# variant calling snake

rule sniffles:
    input:
        initial = "data/assemblies/" + config["assembly"] + ".fasta",
        bam = config["assembly"] + "/outputs/initial/initial_asm.sorted.bam"
    output:
        config["assembly"] + "outputs/variant_calling/" + config["assembly"]
    params:
        threads = config["threads"]
    shell:
        "sniffles -i {~input[1]}  --vcf {output} --reference {input[0]}


# TODO add in some vcf parsing like vcftools or similar. something that can generate stats fairly easily and quickly. 
# doesnt have to be anything amazingm jsut something that can be plotted eventually

# TODO
# one major output of this snake needs to be the alleles plotted against the length of the sacffolds similarly to how i did it for the depth plots
