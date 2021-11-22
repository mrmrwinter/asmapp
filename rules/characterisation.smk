#####GENOME CHARACTERISATION #####

configfile: "config.yaml"

# rule all:
#     input:
#         expand("data/catted/{sample}.{insert}.{FR}P.catted.fq.gz", sample=config["sample"], insert=config["insert"], FR=config["FR"]),
#         expand("reports/kmc/kmer_k21.hist"),
#         expand("reports/genomescope/plot.png"),
#         expand("reports/smudge/smudgeplot_smudgeplot.png")


#concatenation of trimmed reads
rule concatenation:
    conda:
        "../envs/characterisation.yaml"
    input:
        P="data/trimmed/{sample}.{insert}.trimmed_{FR}P.fq.gz",
        U="data/trimmed/{sample}.{insert}.trimmed_{FR}U.fq.gz"
    output:
        "data/catted/{sample}.{insert}.{FR}P.catted.fq.gz"
    threads:
        28
    shell:
        "cp {input.P} {output} && cat {input.U} >> {output}"

### KMC
# counts kmers and coverage#uses all trimmed files of all libraries
rule create_kmcfilelist:
    conda:
        "../envs/characterisation.yaml"
    output:
        "reports/kmc/{sample}/kmcfilelist"
    shell:
        "ls data/raw/{sample} | sed -e 's/^/data\/raw\//' > {output}"


rule kmc_count:
    conda:
        "../envs/characterisation.yaml"
    input:
        "reports/kmc/{sample}/kmcfilelist"
    output:
        pre = "reports/kmc/{sample}/kmer_counts.kmc_pre",
        suf = "reports/kmc/{sample}/kmer_counts.kmc_suf"
    threads:
        28
    params:
        kmc = "reports/kmc/{sample}/kmer_counts"
    shell:
        "rm -rf tmp/ && \
        mkdir tmp/ && \
        kmc \
        -k21 \
        -t{threads} \
        -m50 \
        -ci1 \
        -cs10000 \
        @{input} \
        {params.kmc} \
        tmp/"


rule kmc_transform:
    conda:
        "../envs/characterisation.yaml"
    input:
        suf = "reports/kmc/{sample}/kmer_counts.kmc_suf",
        pre = "reports/kmc/{sample}/kmer_counts.kmc_pre"
    output:
        "reports/kmc/{sample}/kmer_k21.hist"
    params:
        kmc="reports/kmc/{sample}/kmer_counts"
    shell:
        "kmc_tools transform {params.kmc} histogram {output} -cx10000"

#this rule replaces tabs with spaces to allow running from kmc3 into genomescope
rule kmc2genomescope_transformation:
    input:
        "reports/kmc/{sample}/kmer_k21.hist"
    output:
        "reports/kmc/{sample}/kmer_k21.histo"
    shell:
        "expand -t 1 {input} > {output}"

### GenomeScope
#this is to see kmer spectra, estimate genome size, etc
rule genomescope:
    conda:
        "../envs/characterisation.yaml"
    input:
        "reports/kmc/{sample}/kmer_k21.histo"
    output:
        "reports/genomescope/{sample}/plot.log.png"
    params:
        outdir="reports/genomescope/{sample}/"
    shell:
        "Rscript scripts/genomescope.R {input} 21 150 {params.outdir} 1000 1"

#smudgeplot for predicting ploidy
rule smudgeplot:
    conda:
        "../envs/characterisation.yaml"
    input:
        "reports/kmc/{sample}/kmer_k21.hist"
    output:
        "reports/smudge/{sample}/smudgeplot_smudgeplot.png"
    params:
        counts = "reports/kmc/{sample}/kmer_counts",
        dump = "reports/smudge/{sample}/kmer_k21.dump",
        pairs = "reports/smudge/{sample}/kmer_pairs",
        cov = "reports/smudge/{sample}/kmer_pairs_coverages.tsv"
    shell:
        """
        L=$(smudgeplot.py cutoff {input} L)
        U=$(smudgeplot.py cutoff {input} U)

        echo $L $U

        kmc_tools transform {params.counts} \
        -ci$L \
        -cx$U dump \
        -s {params.dump}

        smudgeplot.py hetkmers \
        -o {params.pairs} < {params.dump} \

        smudgeplot.py plot {params.cov}
        """
