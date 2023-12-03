rule count_kmer_in_reads:
    input:
        fq1=getfq1,
        fq2=getfq2
    output: "results/read_kmers/{sample}.kff" 
    params:
        wdir="temp_kmc_{sample}",
        filel="temp_kmc_filelist_{sample}.txt"
    threads: 8
    benchmark: 'benchmark/results/read_kmers/{sample}.count_kmer_in_reads.benchmark.tsv'
    container: "docker://quay.io/biocontainers/kmc:3.2.1--hf1761c0_2"
    shell:
        """
        echo {input.fq1} > {params.filel}
        echo {input.fq2} >> {params.filel}
        mkdir -p {params.wdir}
        kmc -k29 -m64 -okff -t{threads} @{params.filel} {params.wdir}/out {params.wdir}
        mv {params.wdir}/out.kff {output}
        rm -r {params.filel} {params.wdir}
        """

rule map_short_reads_giraffe:
    input: 
        fq1=getfq1,
        fq2=getfq2,
        gbz="results/sample_pg/{graph}.{sample}.gbz",
        dist='results/sample_pg/{graph}.{sample}.dist',
        min='results/sample_pg/{graph}.{sample}.min'
    output: "results/gaf/{sample}.{graph}.gaf.gz"
    threads: 8
    benchmark: 'benchmark/results/gaf/{sample}.{graph}.map_short_reads_giraffe.benchmark.tsv'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    shell:
        """
        vg giraffe --progress \
        --read-group "ID:1 LB:lib1 SM:{wildcards.sample} PL:illumina PU:unit1" \
        --sample "{wildcards.sample}" \
        --output-format gaf \
        -f {input.fq1} -f {input.fq2} \
        -Z {input.gbz} \
        -d {input.dist} \
        -m {input.min} \
        -t {threads} | gzip > {output}
        """

rule sample_haplotypes:
    input: 
        gbz=getgbz(),
        hapl=gethapl(),
        read_kmer='results/read_kmers/{sample}.kff'
    output: "results/sample_pg/{graph}.{sample}.gbz"
    threads: 8
    benchmark: 'benchmark/results/sample_pg/{graph}.{sample}.sample_haplotypes.benchmark.tsv'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    shell:
        """
        vg haplotypes -v 2 -t {threads} \
        --num-haplotypes 4 \
        --present-discount 0.9 \
        --het-adjustment 0.05 \
        --absent-score 0.8 \
        --include-reference \
        --diploid-sampling \
        -i {input.hapl} \
        -k {input.read_kmer} \
        -g {output} {input.gbz}
        """

rule surject_reads:
    input:
        gbz=getgbz(),
        paths_list=config['ref_paths_list'],
        gaf="results/gaf/{sample}.{graph}.gaf.gz"
    output: "results/bam/{sample}.{graph}.bam"
    params:
        half_threads=lambda wildcards, threads: int(threads/2)
    threads: 8
    benchmark: 'benchmark/results/bam/{sample}.{graph}.surject_reads.benchmark.tsv'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    shell:
        """
        vg surject \
        -F {input.paths_list} \
        -x {input.gbz} \
        -t {params.half_threads} \
        --bam-output --gaf-input \
        --sample {wildcards.sample} \
        --read-group "ID:1 LB:lib1 SM:{wildcards.sample} PL:illumina PU:unit1" \
        --prune-low-cplx --interleaved --max-frag-len 3000 \
        {input.gaf} | samtools sort --threads {params.half_threads} \
        -O BAM > {output}
        """
