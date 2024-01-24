# GENOME CHARACTERISATION
# Perform genome profiling with genomescope and smudgeplot

# Count the kmers in the raw reads
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
        threads = config["threads"],
        log = f"{config['assembly']}/logs/{rule}.log",
    shell:
        "rm -rf {params[1]} && \
        mkdir -p {params[1]} && \
        kmc \
        -k21 \
        -t{params[threads]} \
        -m50 \
        -ci1 \
        -cs10000 \
        {input[reads]} \
        {params[kmc]} \
        {params[tmp]} 2> {params[log]}"


# Transform them and generate a histogram of the distribution
rule kmc_transform:
    # conda:
    #     "../envs/characterisation.yaml"
    input:
        suf = config["assembly"] + "/reports/kmc/kmer_counts.kmc_suf",
        pre = config["assembly"] + "/reports/kmc/kmer_counts.kmc_pre"
    output:
        config["assembly"] + "/reports/kmc/kmer_k21.hist"
    params:
        kmc = config["assembly"] + "/reports/kmc/kmer_counts",
        log = f"{config['assembly']}/logs/{rule}.log",
    shell:
        "kmc_tools transform {params[kmc]} histogram {output} -cx10000 2> {params[log]}"


# Replaces tabs with spaces to allow running from kmc3 into genomescope
rule kmc2genomescope_transformation:
    input:
        hist = config["assembly"] + "/reports/kmc/kmer_k21.hist"
    output:
        config["assembly"] + "/reports/kmc/kmer_k21.histo"
    params:
        log = f"{config['assembly']}/logs/{rule}.log",
    shell:
        "expand -t 1 {input[hist]} > {output} 2> {params[log]}"




### GenomeScope
# Run Genomescope to visualise kmer spectra, estimate genome size, etc
rule genomescope:
    # conda:
    #     "../envs/characterisation.yaml"
    input:
        histo = config["assembly"] + "/reports/kmc/kmer_k21.hist"
    output:
        report(
            config["assembly"] + "/reports/genomescope/plot.png",
            caption="../docs/captions/genomescope.rst",
            category="Genome profiling"
        )
    params:
        outdir=config["assembly"] + "/reports/genomescope/",
        ploidy = config["ploidy"],
        log = f"{config['assembly']}/logs/{rule}.log",
    shell:
        "genomescope.R {input} 21 15000 {params.outdir} 1000 1 -p {params[ploidy]} 2> {params[log]}"


# Run smudgeplot to predict ploidy
rule smudgeplot:
    # conda:
    #     "../envs/characterisation.yaml"
    input:
        hist = config["assembly"] + "/reports/kmc/kmer_k21.hist"
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
        cov = config["assembly"] + "/reports/smudge/kmer_pairs_coverages.tsv",
        log = f"{config['assembly']}/logs/{rule}.log",
    shell:
        """
        L=$(smudgeplot.py cutoff {input[hist]} L)
        U=$(smudgeplot.py cutoff {input[hist]} U)

        echo $L $U

        kmc_tools transform {params[counts]} \
        -ci$L \
        -cx$U dump \
        -s {params[dump]}  2> {params[log]}

        smudgeplot.py hetkmers \
        -o {params[pairs]} < {params[dump]} 2>> {params[log]}

        smudgeplot.py plot {params[cov]} 2>> {params[log]}
        """
