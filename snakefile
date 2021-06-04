# KARYON LONGREAD

configfile: "config.yaml"

ruleorder: index_mapping > samtools_index > samtools_faidx > GATK > bcftools > samtools_flagstats

rule all:
    input:
        plot = config["assembly"] + "/outputs/plots/plot.png",
        flagstats = config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.flagstat",
        nucmer = config["assembly"] + "/reports/nucmer.initial.plot.rplot"

rule redundans:
    conda:
        "envs/redundans.yaml"
    input:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",
        reads = "data/reads/" + config["reads"] + ".fastq.gz"
    output:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fa",
    params:
        out_pfx = config["assembly"] + "/outputs/redundans/"
    shell:
        "rm -rf {params[out_pfx]} && " + config["redundans_path"] + "/redundans.py --fasta {input[assembly]} --outdir {params[out_pfx]} --threads 12 --longreads {input[reads]} --verbose"


rule redundans_tagging:
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fa",
    output:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta",
    run:
        from Bio import SeqIO
        to_add = "redundans"
        with open(output[0], "w") as outputs:
            for r in SeqIO.parse(input[0], "fasta"):
                r.id = (to_add + r.description).replace(" ", "_")
                r.description = r.id
                SeqIO.write(r, outputs, "fasta")


rule dictionary_creation:
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta"
    output:
        dictionary = config["assembly"] + "/outputs/redundans/scaffolds.reduced.dict"
    shell:
       "picard CreateSequenceDictionary R={input[assembly]} O={output}"

rule index_mapping:
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta"
    output:
        fai = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fa.fai"
    shell:
       "samtools faidx {input[assembly]} > {output}"        # change this to minimap

rule mapping:
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta",
        reads = "data/reads/" + config["reads"] + ".fastq.gz"
    output:
        sam = config["assembly"] + "/outputs/mapping/scaffolds.reduced.sam"
    params:
        threads = config["threads"],
        seq_tech = "map-ont"
    shell:
        "minimap2 -t {params[threads]} -ax {params[seq_tech]} {input[assembly]} {input[reads]} > {output}"


rule conversion:
    input:
        sam = config["assembly"] + "/outputs/mapping/scaffolds.reduced.sam"
    output:
        bam = config["assembly"] + "/outputs/mapping/scaffolds.reduced.bam"
    shell:
        "samtools view -b -S {input} > {output}"

rule sorting:
    input:
        config["assembly"] + "/outputs/mapping/scaffolds.reduced.bam"
    output:
        bam = config["assembly"] + "/outputs/mapping/scaffolds.reduced.sorted.bam"
    shell:
        "samtools sort {input} > {output}"

rule nucmer_reduced_vs_initial:
    conda:
        "envs/redundans.yaml"
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta",
        initial = "data/assemblies/" + config["assembly"] + ".fasta"
    output:
        config["assembly"] + "/reports/nucmer/nucmer.initial.plot.rplot",
    shell:
        """
        mkdir -p tmp/
        nucmer -p tmp/nucmer.contigs {input[0]} {input[1]}
        mummerplot --png --large tmp/nucmer.contigs.delta -p {output}
        rm -rf tmp/
        """

rule samtools_index:
    input:
       bam = config["assembly"] + "/outputs/mapping/scaffolds.reduced.sorted.bam"
    output:
       bai = config["assembly"] + "/outputs/mapping/scaffolds.reduced.sorted.bam.bai"
    shell:
       "samtools index {input[bam]} > {output}"

rule samtools_faidx:
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta"
    output:
        fai = config["assembly"] + "/outputs/mapping/scaffolds.reduced.fa.fai"
    shell:
        "samtools faidx {input[assembly]} > {output}"

rule GATK:
    container:
        "docker://broadinstitute/gatk:4.0.2.0"
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta",
        bam = config["assembly"] + "/outputs/mapping/scaffolds.reduced.sorted.bam",
        dict = config["assembly"] + "/outputs/redundans/scaffolds.reduced.dict"
    output:
        vcf = config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.vcf"
    params:
        memory = config["memory"] + "G"
    shell:
        "gatk --java-options '-Xmx{params[memory]}' HaplotypeCaller -R {input[assembly]} -I {input[bam]} -O {output}"

rule bcftools:
    conda:
        "envs/v_calling.yaml"
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta",
        bam = config["assembly"] + "/outputs/mapping/scaffolds.reduced.sorted.bam"
    output:
        mpileup = config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.mpileup"
    shell:
        "bcftools mpileup --fasta-ref {input[assembly]} {input[bam]} > {output}"

rule samtools_flagstats:
    input:
        bam = config["assembly"] + "/outputs/mapping/scaffolds.reduced.sorted.bam"
    output:
        flagstats = config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.flagstat"
    shell:
        "samtools flagstats {input[bam]} > {output}"

rule karyon_plots:
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta",
        mpileup = config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.mpileup",
        flagstats = config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.flagstat",
        bam = config["assembly"] + "/outputs/mapping/scaffolds.reduced.sorted.bam",
        vcf = config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.vcf"
    output:
        plot = config["assembly"] + "/outputs/plots/plot.png"
    params:
        out_pfx = config["assembly"] + "/outputs/plots",
        out_name = config["assembly"]
    shell:
        "python3 scripts/karyon_plots.py --fasta {input[assembly]} --output_directory {params[out_pfx]} --output_name {params[out_name]} --vcf {input[vcf]} --pileup {input[mpileup]} --bam {input[bam]} --library --configuration --wsize --max_scaf2plot --scafminsize --scafmaxsize --job_id"        # Identifier of the intermediate files generated by the different programs. If false, the program will assign a name consisting of a string of 6 random alphanumeric characters.')"


------------------------------

rule make_blast_database:  # Rule to make database of cds fasta
    input:
        config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta" # input to the rule
    output:
        nhr = "data/database/" + config["assembly"] + ".nhr",   # all outputs expected from the rule
        nin = "data/database/" + config["assembly"] + ".nin",
        nsq = "data/database/" + config["assembly"] + ".nsq"
    params:
        "data/database/" + config["assembly"]   # prefix for the outputs, required by the command
    shell:  # shell command for the rule
        "makeblastdb \
        -in {input} \
        -out {params} \
        -dbtype nucl"  # the database type


rule blast:
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta",
        nin = "data/database/" + config["assembly"] + ".nin"
    output:
        tsv
    run:
        from Bio import SeqIO
        import pandas as pd

        seqs = input[0]
        tsv = output[0]

        for contig in seqs:
            blast against assembly database
            Take top hit row that is not self
            put the entry for that hit into a pandas frame

        print pandas frame to tsv


rule dotplots:
    input:
        tsv,
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta"
    output:
        dotplots
    run:
        open tsv as dataframe
        open assembly with bio

        for row in dataframe:
            q = column1 seq pulled with Bio
            h = column2 seq pulled with Bio
            Dgenies q against h
