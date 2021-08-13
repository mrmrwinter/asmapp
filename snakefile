# assembly appraisal (karyinon workflow)


configfile: "config.yaml"


###############################################################################

rule all:
    input:
# karyon
        flagstats = config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.flagstat",
        vcf = config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.vcf",
        mpileup = config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.mpileup",
        plot = config["assembly"] + "/outputs/plots/plot.png",
# dotplots
    # NUCMER
        nucmer = config["assembly"] + "/reports/nucmer/nucmer.initial.png",
        nucmer_ref = config["assembly"] + "/reports/nucmer/nucmer.reference.png",
        nucmer_ref_int = config["assembly"] + "/reports/nucmer/nucmer.int_ref.png",
    # BLAST
        tsv = config["assembly"] + "/reports/blast/blast.out",
# assembly stats
#        quast = config["assembly"] + "/reports/quast/report.txt",



###############################################################################

# ASSEMBLY COLLAPSING

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





###############################################################################

# READ MAPPING AND MANAGEMENT

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



rule mapping:
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta",
        reads = "data/reads/" + config["reads"] + ".fastq.gz"
    output:
        sam = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sam"
    params:
        threads = config["threads"],
        seq_tech = "map-" + config["seq_tech"]
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


rule samtools_index:
    input:
       bam = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sorted.bam"
    output:
       bai = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sorted.bam.bai",
    shell:
       "samtools index {input[bam]} > {output[0]}"


rule samtools_faidx:
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta"
    output:
        fai = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta.fai"
    shell:
        "samtools faidx {input[assembly]}"


rule samtools_tagged_index:
    input:
       bam = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sorted.tagged.bam"
    output:
       fai = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sorted.tagged.bam.bai",
       nuc_fai = config["assembly"] + "/reports/nucmer/scaffolds.reduced.fasta.fai"
    shell:
       "samtools index {input[bam]} && cp {output[0]} {output[1]}"



###############################################################################

# NUCMER ASSEMBLY DOTPLOTS

rule nucmer_reduced_vs_initial:
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta",
        initial = "data/assemblies/" + config["assembly"] + ".fasta",
        nuc_fai = config["assembly"] + "/reports/nucmer/scaffolds.reduced.fasta.fai"
    output:
        config["assembly"] + "/reports/nucmer/nucmer.initial.png",
        config["assembly"] + "/reports/nucmer/nucmer.initial.delta"
    params:
        "nucmer.initial",
        config["assembly"] + "/reports/nucmer/",
    shell:
        """
        mkdir -p tmp/
        nucmer -p tmp/{params[0]} {input[0]} {input[1]}
        cp tmp/{params[0]}.delta {params[1]}
        mummerplot -l -f --png --large {params[1]}{params[0]}.delta -p {params[1]}{params[0]}
        """


rule nucmer_reduced_vs_reference:
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta",
        reference = config["reference"] + ".fasta.gz",
        nuc_fai = config["assembly"] + "/reports/nucmer/scaffolds.reduced.fasta.fai"
    output:
        config["assembly"] + "/reports/nucmer/nucmer.reference.png",
        config["assembly"] + "/reports/nucmer/nucmer.reference.delta"
    params:
        "nucmer.reference",
        config["assembly"] + "/reports/nucmer/",
        config["reference"] + ".fasta"
    shell:
        """
        mkdir -p tmp/
        gunzip -c {input[reference]} > {params[2]}
        nucmer -p tmp/{params[0]} {input[0]} {params[2]}
        cp tmp/{params[0]}.delta {params[1]}
        mummerplot -l -f --png --large {params[1]}{params[0]}.delta -p {params[1]}{params[0]}
        """


rule nucmer_initial_vs_reference:
    input:
        initial = "data/assemblies/" + config["assembly"] + ".fasta",
        reference = config["reference"] + ".fasta.gz",
        nuc_fai = config["assembly"] + "/reports/nucmer/scaffolds.reduced.fasta.fai"
    output:
        config["assembly"] + "/reports/nucmer/nucmer.int_ref.png",
        config["assembly"] + "/reports/nucmer/nucmer.int_ref.delta"
    params:
        "nucmer.int_ref",
        config["assembly"] + "/reports/nucmer/",
        config["reference"] + ".fasta"
    shell:
        """
        mkdir -p tmp/
        gunzip -c {input[reference]} > {params[2]}
        nucmer -p tmp/{params[0]} {input[0]} {params[2]}
        cp tmp/{params[0]}.delta {params[1]}
        mummerplot -l -f --png --large {params[1]}{params[0]}.delta -p {params[1]}{params[0]}
        """


 # maybew better to add the circle plots as an addition to the other nucmer rules to prevent having 6 rules
 # know for sure that the nucmer rules currently work
 #  can tag on shell command to execute the script
  # have nucmer delta as an output then pull that output into the rscript command to generate the second output, the circle plots

 # set up the outputs to generate the delta files for the circle plots but i
 # think i need to keep it as a script execuition in order for the snakemake inputs to have access to it


# rule nucmer_circles:
#     input:
#         fai = config["assembly"] + "/reports/nucmer/scaffolds.reduced.fasta.fai"
#     output:
#         config["assembly"] + "/reports/nucmer/circle.int_ref.png",
#     script:
#         "scripts/snmk_nucmer_circles.R"





###############################################################################

# KARYON PLOIDY PLOTS

rule dictionary_creation:
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta"
    output:
        dictionary = config["assembly"] + "/outputs/redundans/scaffolds.reduced.dict"
    shell:
       "picard CreateSequenceDictionary R={input[assembly]} O={output}"


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
#    conda:
#        "envs/v_calling.yaml"
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta",
        faidx = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta.fai",
        bam = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sorted.bam",
        bai = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sorted.bam.bai",
    output:
        mpileup = config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.mpileup"
    shell:
        # "bcftools mpileup -Ou -f {input[0]} {input[2]} | \
        # bcftools call -Ou -mv | \
        # bcftools filter -s LowQual -e '%QUAL<10 || DP>100' > {output}"
        # "bcftools mpileup -Ov -f {input[assembly]} {input[bam]} -o {output}"
        "samtools mpileup -B -f {input[assembly]} -v {input[bam]} -o {output}"

rule samtools_flagstats:
    input:
        bam = config["assembly"] + "/outputs/redundans/scaffolds.reduced.sorted.tagged.bam"
    output:
        flagstats = config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.flagstat"
    shell:
        "samtools flagstats {input[bam]} > {output}"


rule karyon_plots:
    conda:
        "envs/karyonplot.yaml"
    input:
        assembly = config["full_path"] + config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta",
        mpileup = config["full_path"] + config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.mpileup",
        flagstats = config["full_path"] + config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.flagstat",
        bam = config["full_path"] + config["assembly"] + "/outputs/redundans/scaffolds.reduced.sorted.tagged.bam",
        vcf = config["full_path"] + config["assembly"] + "/outputs/variant_calling/scaffolds.reduced.vcf",
        reads = config["full_path"] + "data/reads/" + config["reads"] + ".fastq.gz"
    output:
        plot = config["assembly"] + "/outputs/plots/plot.png"
    params:
        out_pfx = config["assembly"] + "/outputs/plots",
        out_name = config["assembly"]
    shell:
        "python3 scripts/karyon_plots.py --fasta {input[assembly]} --output_directory {params[out_pfx]} --output_name {params[out_name]} --vcf {input[vcf]} --pileup {input[mpileup]} --bam {input[bam]} --library {input[reads]}" # --configuration -wsize --max_scaf2plot --scafminsize --scafmaxsize --job_id"        # Identifier of the intermediate files generated by the different programs. If false, the program will assign a name consisting of a string of 6 random alphanumeric characters.')"



##############################################################################

# PAIRS TABLE with BLAST

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
        db = "data/database/" + config["assembly"] + "/" + config["assembly"] + ".nin"
    output:
        tsv = config["assembly"] + "/reports/blast/blast.out"
    params:
        out_pfx = config["assembly"] + "/reports/blast",
        db_pfx = "data/database/" + config["assembly"] + "/" + config["assembly"]
    shell:
        "blastn -query {input[0]} -db {params[1]} -outfmt 6 -max_target_seqs 2 -out {params[0]}/blast.out"


##############################################################################

# PAIRS TABLE WITH LASTZ

rule lastz:
    input:
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta"
    output:
        tsv = config["assembly"] + "/reports/lastz/
    params:
        out_pfx = config["assembly"] + "/reports/lastz"
    shell:




###############################################################################

# QUAST

rule reads_to_fasta:
    input:
        reads = "data/reads/" + config["reads"] + ".fastq.gz",
    output:
        reads = "data/reads/" + config["reads"] + ".fasta",
    shell:
        "zcat -c {input} | seqkit fq2fa | cat > {output}"


# rule download_busco_for_quast:
#     message:
#         "[INFO] Downloading BUSCO databases for QUAST appraisal..."
#     output:
#         ""
#     shell:
#         "quast-download-busco"


  # Performing QUAST assembly appraisal
rule quast:
    message:
        "[INFO] Performing QUAST appraisal on assemblies..."
    input:
        initial = "data/assemblies/" + config["assembly"] + ".fasta",
        assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fa",
        reads = "data/reads/" + config["reads"] + ".fasta",
        reference = config["reference"] + ".fasta.gz"
    output:
        config["assembly"] + "/reports/quast/report.txt",
    params:
        out_pfx = config["assembly"] + "/reports/quast/",
        threads = config["threads"]
    shell:
        config["quast_path"] + "/quast.py --large {input[0]} {input[1]} --glimmer -b --threads {params[1]} -L -r {input[2]} --nanopore {input[reads]} -o {params[0]}"


##############################################################################

# BLOBPLOTS

# rule blobs:
#     input:
#         assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fa",
#     output:
#     shell:
#         """
#         """
#
#



##############################################################################

# CEGMA

# rule CEGMA_collapsed:
#     container:
#         "chrishah/cegma:2.5"
#     input:
#         assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fa",
#     output:
#     shell:
#
#



##############################################################################

# MITOCHONDRIAL DETECTION

# rule mito_identification:
#     input:
#         initial = "data/assemblies/" + config["assembly"] + ".fasta",
#
#



##############################################################################

# COVERAGE PLOTS

# rule coverage_plots:
#     input:
#         initial = "data/assemblies/" + config["assembly"] + ".fasta",




##############################################################################

# MERQURY

    # Performing Merqury assembly appraisal
# rule merqury:
#    message:
#        "[INFO] Performing Merqury assembly appraisal..."
#    conda:
#        "../envs/merqury.yaml"
#    input:
#        initial = "data/assemblies/" + config["assembly"] + ".fasta",
#    output:
#    params:
#    shell:













###############################################################################
# FLUFF
##################################################

# rule dotplots:
#     input:
#         tsv = config["assembly"] + "/reports/mapped_nonself_hits.sam",
#         assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta"
#     output:
#         config["assembly"] + "/reports/dotplots/dotplot.png"
#     params:
#         out_dir = config["assembly"] + "/reports/dotplots/"
#     run:
#         from Bio import SeqIO
#         import pandas as pd
#         import os
#         import sys
#         from readpaf import parse_paf
#         import dotplot
#
#         shell("mkdir -p tmp/")
#
#         tmp_tbl = input[0].replace(".paf",".cleaned.sam")
#         print(tmp_tbl)
#         shell("sed '/^@/ d' < " + input[0] + " > " + tmp_tbl)
#
#         df = pd.read_csv(tmp_tbl, sep = '\t')
#         seqs = input[1]
# #
#         for index, value in df.iterrows():
#             query = value[0]
#             hit = value[2]
#             print(query + " pairs with " + hit)
#             for record in SeqIO.parse(seqs,'fasta'):
#                 q_seq=[]
#                 h_seq=[]
#                 if record.id == query:
#                     q_seq = ">" + record.id + "\n" + record.seq + "\n"
#                     with open('tmp/' + query + '.fasta', 'a+', newline='\n') as g:
#                         SeqIO.write(record, g, 'fasta-2line')
#                     g.close()
#                     print(q_seq)
#                 if record.id == hit:
#                     h_seq = ">" + record.id + "\n" + record.seq + "\n"
#                     with open('tmp/' + hit + '.fasta', 'a+', newline='\n') as r:
#                         SeqIO.write(record, r, 'fasta-2line')
#                     r.close()
#                     print(h_seq)
#             shell("dotplot --drawer matplotlib --fasta tmp/" + str(query) + ".fasta tmp/" + str(hit) + ".fasta > " + params[0] + str(query) + ".dotplot.png")



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

# rule mapping_back:
#     input:
#         assembly = config["assembly"] + "/outputs/redundans/scaffolds.reduced.fasta"
#     output:
#         tsv = config["assembly"] + "/reports/minimap2/mapped_nonself_hits.sam"
#     params:
#         threads = config["threads"]
#     shell:
#         "minimap2 -ax asm5 -D {input[0]} {input[0]} > {output}"
#
# # rule removing_nonself:
# #     input:
# #         tsv = config["assembly"] + "/reports/mapped_hits.sam"
# #     output:
# #         tsv = config["assembly"] + "/reports/mapped_nonself_hits.sam"
# #     shell:
# #         "samtools view -F0x900 {input} > {output}"
