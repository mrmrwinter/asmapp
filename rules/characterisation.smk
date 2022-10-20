# GENOME CHARACTERISATION

### KMC
# counts kmers and coverage
# rule create_kmcfilelist:
#     conda:
#         "../envs/characterisation.yaml"
#     input:
#         reads = "data/reads/" + config["reads"] + ".fastq.gz",
#     output:
#         config["assembly"] + "/reports/kmc/kmcfilelist"
#     shell:
#         "ls {input} | sed -e 's/^data\/reads\///' > {output}"


rule kmc_count:
    # conda:
    #     "../envs/characterisation.yaml"
    input:
        reads = "data/reads/" + config["reads"] + ".fastq.gz",
    output:
        pre = config["assembly"] + "/reports/kmc/kmer_counts.kmc_pre",
        suf = config["assembly"] + "/reports/kmc/kmer_counts.kmc_suf"
    threads:
        28
    params:
        kmc = config["assembly"] + "/reports/kmc/kmer_counts",
        tmp = config["assembly"] + "/reports/kmc/tmp",
        threads = config["threads"]
    shell:
        "rm -rf {params[1]} && \
        mkdir -p {params[1]} && \
        kmc \
        -k21 \
        -t{params[2]} \
        -m50 \
        -ci1 \
        -cs10000 \
        {input} \
        {params[0]} \
        {params[1]}"


rule kmc_transform:
    # conda:
    #     "../envs/characterisation.yaml"
    input:
        suf = config["assembly"] + "/reports/kmc/kmer_counts.kmc_suf",
        pre = config["assembly"] + "/reports/kmc/kmer_counts.kmc_pre"
    output:
        config["assembly"] + "/reports/kmc/kmer_k21.hist"
    params:
        kmc=config["assembly"] + "/reports/kmc/kmer_counts"
    shell:
        "kmc_tools transform {params.kmc} histogram {output} -cx10000"

#this rule replaces tabs with spaces to allow running from kmc3 into genomescope
rule kmc2genomescope_transformation:
    input:
        config["assembly"] + "/reports/kmc/kmer_k21.hist"
    output:
        config["assembly"] + "/reports/kmc/kmer_k21.histo"
    shell:
        "expand -t 1 {input} > {output}"

### GenomeScope
#this is to see kmer spectra, estimate genome size, etc
rule genomescope:
    # conda:
    #     "../envs/characterisation.yaml"
    input:
        config["assembly"] + "/reports/kmc/kmer_k21.histo"
    output:
        report(
            config["assembly"] + "/reports/genomescope/plot.png",
            caption="../docs/captions/genomescope.rst",
            category="Genome profiling"
        )
    params:
        outdir=config["assembly"] + "/reports/genomescope/"
    shell:
        "Rscript scripts/genomescope.R {input} 21 150 {params.outdir} 1000 1"


#smudgeplot for predicting ploidy
rule smudgeplot:
    # conda:
    #     "../envs/characterisation.yaml"
    input:
        config["assembly"] + "/reports/kmc/kmer_k21.hist"
    output:
        report(
            config["assembly"] + "/reports/smudge/smudgeplot_smudgeplot.png",
            caption="../docs/captions/smudgeplot.rst",
            category="Genome profiling"
        )
    params:
        counts = config["assembly"] + "/reports/kmc/kmer_counts",
        dump = config["assembly"] + "/reports/smudge/kmer_k21.dump",
        pairs = config["assembly"] + "/reports/smudge/kmer_pairs",
        cov = config["assembly"] + "/reports/smudge/kmer_pairs_coverages.tsv"
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
