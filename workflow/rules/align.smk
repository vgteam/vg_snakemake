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
        surj_threads=lambda wildcards, threads: max(1, int(threads/2)) if threads < 5 else threads - 2,
        sort_threads=lambda wildcards, threads: max(1, int(threads/2)) if threads < 5 else 2,
        sort_dir="{temp_bam_sort_{sample}_{graph}"
    threads: 8
    benchmark: 'benchmark/results/bam/{sample}.{graph}.surject_reads.benchmark.tsv'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    shell:
        """
        mkdir -p {params.sort_dir}

        vg surject \
        -F {input.paths_list} \
        -x {input.gbz} \
        -t {params.surj_threads} \
        --bam-output --gaf-input \
        --sample {wildcards.sample} \
        --read-group "ID:1 LB:lib1 SM:{wildcards.sample} PL:illumina PU:unit1" \
        --prune-low-cplx --interleaved --max-frag-len 3000 \
        {input.gaf} | samtools sort --threads {params.sort_threads} -T {params.sort_dir}/temp \
        -O BAM > {output}

        rm -rf {params.sort_dir}
        """

rule extract_unmapped_reads:
    input:
        fq1=getfq1,
        fq2=getfq2,
        gaf="results/gaf/{sample}.{graph}.gaf.gz"
    output:
        fq1="results/unmapped_reads/{sample}.{graph}.unmapped.1.fq.gz",
        fq2="results/unmapped_reads/{sample}.{graph}.unmapped.2.fq.gz"
    params:
        reads="temp.{sample}.{graph}.unmapped.reads.txt"
    threads: 1
    benchmark: 'benchmark/results/unmapped_reads/{sample}.{graph}.extract_unmapped_reads.benchmark.tsv'    
    shell:
        """
        zcat {input.gaf} | awk '{{if($3=="*"){{print $1}}}}' | uniq > {params.reads}
        seqtk subseq {input.fq1} {params.reads} | gzip > {output.fq1}
        seqtk subseq {input.fq2} {params.reads} | gzip > {output.fq1}
        rm {params.reads}
        """

rule trim_fastq:
    input:
        fq1=getfq1,
        fq2=getfq2,
        refsynt=config['refsynt_fa'],
        adpt=config['adapters_fa']
    output:
        fq12="results/fq/{sample}.trimmed.fq.gz",
        stats_synt="results/fq/qc.{sample}/bbmap.{sample}.statsSYNT.txt",
        synt_fq1="results/fq/qc.{sample}/{sample}.synt.1.fq.gz",
        synt_fq2="results/fq/qc.{sample}/{sample}.synt.2.fq.gz",
        qc_before="results/fq/qc.{sample}/",
        qc_after="results/fq/qc.{sample}/"
    params:
        qcdir="results/fq/fastqc.{sample}"
    resources:
        mem="20G",
        runtime="3h"
    threads: 4
    shell:
        """
        mkdir {params.qcdir}
        fastqc --noextract -t {threads} -o {params.qcdir} -a {input.adapters_fa} -c {input.refsynt_fa} {input.fq1} {input.fq2}
        
        seqtk mergepe {input.fq1} in2={input.fq2} | \
        bbduk.sh -Xmx10g -t={threads} in=stdin.fq out=stdout.fq interleaved=t ref={input.adapters_fa} \
        ftm=5 minlen=25 qtrim=rl trimq=10 ktrim=r k=23 mink=11 hdist=1 tpe tbo | \
        bbduk.sh -Xmx10g -t={threads} in=stdin.fq out={output.fq12} interleaved=t \
        outm1={output.synt_fq1} outm2={output.synt_fq1} ref={input.refsynt_fa} k=31 hdist=1 stats={output.stats_synt}

        fastqc --noextract -t {threads} -o {params.qcdir} -a {input.adapters_fa} -c {input.refsynt_fa} {output.fq12}
        """
