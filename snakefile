# KARYON LONGREAD

configfile: "config.yaml"
#
# ruleorder: index_mapping > samtools_index > samtools_faidx > GATK > bcftools > samtools_flagstats

rule all:
    input:
        plot = config["assembly"] + "/outputs/plots/plot.png",
        flagstats = config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.flagstat",
        nucmer = config["assembly"] + "/reports/nucmer/nucmer.initial.rplot",
        ex_tsv = config["assembly"] + "/reports/minimap2/mapped_nonself_hits.sam",
        # out_dir = config["assembly"] + "/reports/dotplots/dotplot.png",
        tsv = config["assembly"] + "/reports/blast/blast.out",
        quast = config["assembly"] + "/reports/quast/report.txt",


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
        to_add = "redundans_"
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

# rule index_mapping:
#     input:
#         assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta"
#     output:
#         fai = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fa.fai"
#     shell:
#        "samtools faidx {input[assembly]}"        # change this to minimap

rule mapping:
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta",
        reads = "data/reads/" + config["reads"] + ".fastq.gz"
    output:
        sam = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sam"
    params:
        threads = config["threads"],
        seq_tech = "map-ont"
    shell:
        "minimap2 -t {params[threads]} -ax {params[seq_tech]} {input[assembly]} {input[reads]} > {output}"


rule conversion:
    input:
        sam = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sam"
    output:
        bam = config["assembly"] + "/outputs/redundans/scaffolds.reduced.bam"
    shell:
        "samtools view -b -S {input} > {output}"

rule sorting:
    input:
        config["assembly"] + "/outputs/redundans/scaffolds.reduced.bam"
    output:
        bam = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sorted.bam"
    shell:
        "samtools sort {input} > {output}"

rule nucmer_reduced_vs_initial:
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta",
        initial = "data/assemblies/" + config["assembly"] + ".fasta"
    output:
        config["assembly"] + "/reports/nucmer/nucmer.initial.rplot",
    params:
        config["assembly"] + "/reports/nucmer/nucmer.initial",
        config["assembly"] + "/reports/nucmer/",
    shell:
        """
        mkdir -p tmp/
        nucmer -p tmp/nucmer.contigs {input[0]} {input[1]}
        cp tmp/nucmer.contigs.delta {params[1]}
        mummerplot --png --large {params[1]}nucmer.contigs.delta -p {params[0]}
        rm -rf tmp/
        """

rule samtools_index:
    input:
       bam = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sorted.tagged.bam"
    output:
       bai = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sorted.tagged.bam.bai"
    shell:
       "samtools index {input[bam]}"

rule samtools_faidx:
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta"
    output:
        fai = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta.fai"
    shell:
        "samtools faidx {input[assembly]}"

rule fix_bam:
    input:
        bam = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sorted.bam",
    output:
        bam = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sorted.tagged.bam",
    params:
        seq_tech = config["seq_tech"]
    shell:
        "picard AddOrReplaceReadGroups \
        I={input} \
        O={output} \
        RGLB=lib1 \
        RGPL={params[seq_tech]} \
        RGPU=unit1 \
        RGSM=Random"

rule GATK:
    container:
        "docker://broadinstitute/gatk:4.0.2.0"
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta",
        faidx = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta.fai",
        bam = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sorted.tagged.bam",
        bai = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sorted.tagged.bam.bai",
        dict = config["assembly"] + "/outputs/redundans/scaffolds.reduced.dict"
    output:
        vcf = config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.vcf"
    params:
        memory = config["memory"] + "G",
    shell:
        "gatk --java-options -Xmx{params} HaplotypeCaller -R {input[assembly]} -I {input[bam]} -O {output}"

rule bcftools:
    conda:
        "envs/v_calling.yaml"
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta",
        faidx = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta.fai",
        bam = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sorted.tagged.bam",
        bai = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sorted.tagged.bam.bai",
    output:
        mpileup = config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.mpileup"
    shell:
        # "bcftools mpileup -Ov -f {input[assembly]} {input[bam]} -o {output}"
        "samtools mpileup -f {input[assembly]} {input[bam]} | bcftools view -Ov - > {output}"

rule samtools_flagstats:
    input:
        bam = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sorted.bam"
    output:
        flagstats = config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.flagstat"
    shell:
        "samtools flagstats {input[bam]} > {output}"

rule karyon_plots:
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta",
        mpileup = config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.mpileup",
        flagstats = config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.flagstat",
        bam = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sorted.tagged.bam",
        vcf = config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.vcf"
    output:
        plot = config["assembly"] + "/outputs/plots/plot.png"
    params:
        out_pfx = config["assembly"] + "/outputs/plots",
        out_name = config["assembly"]
    shell:
        "python3 scripts/karyon_plots.py --fasta {input[assembly]} --output_directory {params[out_pfx]} --output_name {params[out_name]} --vcf {input[vcf]} --pileup {input[mpileup]} --bam {input[bam]} --library --configuration --wsize --max_scaf2plot --scafminsize --scafmaxsize --job_id"        # Identifier of the intermediate files generated by the different programs. If false, the program will assign a name consisting of a string of 6 random alphanumeric characters.')"


rule mapping_back:
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta"
    output:
        tsv = config["assembly"] + "/reports/minimap2/mapped_nonself_hits.sam"
    params:
        threads = config["threads"]
    shell:
        "minimap2 -ax asm5 -D {input[0]} {input[0]} > {output}"

# rule removing_nonself:
#     input:
#         tsv = config["assembly"] + "/reports/mapped_hits.sam"
#     output:
#         tsv = config["assembly"] + "/reports/mapped_nonself_hits.sam"
#     shell:
#         "samtools view -F0x900 {input} > {output}"

rule make_blast_database:  # Rule to make database of cds fasta
    input:
        config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta" # input to the rule
    output:
        nhr = "data/database/" + config["assembly"] + "/" + config["assembly"] + ".nhr",   # all outputs expected from the rule
        nin = "data/database/" + config["assembly"] + "/" + config["assembly"] + ".nin",
        nsq = "data/database/" + config["assembly"] + "/" + config["assembly"] + ".nsq"
    params:
        "data/database/" + config["assembly"] + "/" + config["assembly"]   # prefix for the outputs, required by the command
    shell:  # shell command for the rule
        "makeblastdb \
        -in {input} \
        -out {params} \
        -dbtype nucl"  # the database type

rule blast_nonself:
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta",
        db = "data/database/" + config["assembly"] + "/" + config["assembly"] + ".nin",

    output:
        tsv = config["assembly"] + "/reports/blast/blast.out"
    params:
        out_pfx = config["assembly"] + "/reports/blast",
        db_pfx = "data/database/" + config["assembly"] + "/" + config["assembly"]
    shell:
        "blastn -query {input[0]} -db {params[1]} -outfmt 6 -max_target_seqs 2 -out {params[0]}/blast.out"






# rule download_busco_for_quast:
#     message:
#         "[INFO] Downloading BUSCO databases for QUAST appraisal..."
#     output:
#         ""
#     shell:
#         "quast-download-busco"






rule reads_to_fasta:
    input:
        reads = "data/reads/" + config["reads"] + ".fastq.gz",
    output:
        reads = "data/reads/" + config["reads"] + ".fasta",
    shell:
        "zcat -c {input} | seqkit fq2fa | cat > {output}"
  # Performing QUAST assembly appraisal
rule quast:
    message:
        "[INFO] Performing QUAST appraisal on assemblies..."
    input:
        redundans = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta",
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fa",
        reads = "data/reads/" + config["reads"] + ".fasta",
        reference = "/media/mike/WD_4TB/javanica_assemblies/Blanc_Mathieu_2017_mJavanica.fasta.gz"
    output:
        config["assembly"] + "/reports/quast/report.txt",
    params:
        out_pfx = config["assembly"] + "/reports/quast/",
        threads = config["threads"]
    shell:
        config["quast_path"] + "/quast.py {input[0]} {input[1]} --large --glimmer -b --threads {params[1]} -L -r {input[2]} --nanopore {input[reads]} -o {params[0]}"

#     # Performing Merqury assembly appraisal
# rule merqury:
#     message:
#         "[INFO] Performing Merqury assembly appraisal..."
#     conda:
#         "../envs/merqury.yaml"
#     input:
#         assembly = "data/assemblies/{assembler}/{sample}/scaffolds.fasta"
#
#     output:
#     params:
#     shell:







    # run:
    #     from Bio import SeqIO
    #     import os
    #
    #     if sys.version_info[0] < 3:
    #         from StringIO import StringIO
    #     else:
    #         from io import StringIO
    #
    #     for seq_record in SeqIO.parse(input[0], "fasta"):
    #         contig = '>' + seq_record.description + '\n' + str(seq_record.seq)
    # -negative_seqidlist exclude_me
# "blastn -query {input[0]} -db {params[1]} -outfmt '6 qseqid sseqid pident' -out {params[0]}/blast.out"
#                 # q =
#                 # h = df.[1]
#                 # def get_first_pd(condition, df):
#                 #     return df[condition(df)].iloc[0]
#                 # q = str(df.iloc[1]).split('\\t')[0]
#                 # q = q.replace('redundansscaffold','')
#                 # q = int(q)
#                 # h = str(df.iloc[1]).split('\\t')[0]
#                 # h = h.replace('redundansscaffold','')
#                 # h = int(h)
#                 # # hit_row = get_first_pd(lambda x: x(contig != h), df)
#                 # print(h)
#                 # writer.writerow(hit_row)
#
#                 shell("rm " + seq_record.id + "*.out")
#
#         g.close()

#
#
rule dotplots:
    input:
        tsv = config["assembly"] + "/reports/mapped_nonself_hits.sam",
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta"
    output:
        config["assembly"] + "/reports/dotplots/dotplot.png"
    params:
        out_dir = config["assembly"] + "/reports/dotplots/"
    run:
        from Bio import SeqIO
        import pandas as pd
        import os
        import sys
        from readpaf import parse_paf
        import dotplot

        shell("mkdir -p tmp/")

        tmp_tbl = input[0].replace(".paf",".cleaned.sam")
        print(tmp_tbl)
        shell("sed '/^@/ d' < " + input[0] + " > " + tmp_tbl)

        df = pd.read_csv(tmp_tbl, sep = '\t')
        seqs = input[1]
#
        for index, value in df.iterrows():
            query = value[0]
            hit = value[2]
            print(query + " pairs with " + hit)
            for record in SeqIO.parse(seqs,'fasta'):
                q_seq=[]
                h_seq=[]
                if record.id == query:
                    q_seq = ">" + record.id + "\n" + record.seq + "\n"
                    with open('tmp/' + query + '.fasta', 'a+', newline='\n') as g:
                        SeqIO.write(record, g, 'fasta-2line')
                    g.close()
                    print(q_seq)
                if record.id == hit:
                    h_seq = ">" + record.id + "\n" + record.seq + "\n"
                    with open('tmp/' + hit + '.fasta', 'a+', newline='\n') as r:
                        SeqIO.write(record, r, 'fasta-2line')
                    r.close()
                    print(h_seq)
            shell("dotplot --drawer matplotlib --fasta tmp/" + str(query) + ".fasta tmp/" + str(hit) + ".fasta > " + params[0] + str(query) + ".dotplot.png")
