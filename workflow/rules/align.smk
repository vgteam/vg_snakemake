rule count_kmer_in_reads:
    input:
        fq1=getfq1,
        fq2=getfq2
    output: "results/{sample}/{sample}.kff" 
    params:
        wdir="temp_kmc_{sample}",
        filel="temp_kmc_filelist_{sample}.txt"
    threads: 8
    benchmark: 'benchmark/{sample}.count_kmer_in_reads.benchmark.tsv'
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
        gbz="results/{sample}/{graph}.sample_pg.{sample}.gbz",
        dist='results/{sample}/{graph}.sample_pg.{sample}.dist',
        min='results/{sample}/{graph}.sample_pg.{sample}.min'
    output: "results/{sample}/{sample}.{graph}.gaf.gz"
    threads: 8
    benchmark: 'benchmark/{sample}.{graph}.map_short_reads_giraffe.benchmark.tsv'
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
        read_kmer='results/{sample}/{sample}.kff'
    output: "results/{sample}/{graph}.sample_pg.{sample}.gbz"
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
        gbz="results/{sample}/{graph}.sample_pg.{sample}.gbz",
        paths_list=config['ref_paths_list'],
        gaf="results/{sample}/{sample}.{graph}.gaf.gz"
    output: "results/{sample}/{sample}.{graph}.bam"
    params:
        surj_threads=lambda wildcards, threads: max(1, int(threads/2)) if threads < 5 else threads - 2,
        sort_threads=lambda wildcards, threads: max(1, int(threads/2)) if threads < 5 else 2,
        sort_dir="{temp_bam_sort_{sample}_{graph}"
    threads: 8
    benchmark: 'benchmark/{sample}.{graph}.surject_reads.benchmark.tsv'
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
        gaf="results/{sample}/{sample}.{graph}.gaf.gz"
    output:
        fq1="results/{sample}/{sample}.{graph}.unmapped.1.fq.gz",
        fq2="results/{sample}/{sample}.{graph}.unmapped.2.fq.gz"
    params:
        reads="temp.{sample}.{graph}.unmapped.reads.txt"
    threads: 1
    benchmark: 'benchmark/{sample}.{graph}.extract_unmapped_reads.benchmark.tsv'    
    shell:
        """
        zcat {input.gaf} | awk '{{if($3=="*"){{print $1}}}}' | uniq > {params.reads}
        seqtk subseq {input.fq1} {params.reads} | gzip > {output.fq1}
        seqtk subseq {input.fq2} {params.reads} | gzip > {output.fq2}
        rm {params.reads}
        """

# rule trim_fastq:
#     input:
#         fq1=getfq1,
#         fq2=getfq2,
#         refsynt_fa=config['refsynt_fa'],
#         adpt_fa=config['adapters_fa'],
#         adpt_tsv=config['adapters_tsv']
#     output:
#         fq12="results/{sample}/{sample}.trimmed.fq.gz",
#         qc_zip="results/{sample}/fastq_qc.{sample}.zip"
#     params:
#         qcdir="temp_fastqc.{sample}",
#         stats_synt="temp_fastqc.{sample}/{sample}.statsSYNT.txt",
#         synt_fq1="temp_fastqc.{sample}/{sample}.synt.1.fq.gz",
#         synt_fq2="temp_fastqc.{sample}/{sample}.synt.2.fq.gz"
#     benchmark: 'benchmark/{sample}.trim_fastq.benchmark.tsv'    
#     resources:
#         mem="8G",
#         runtime="3h"
#     threads: 4
#     shell:
#         """
#         mkdir {params.qcdir}
#         fastqc --extract --delete -f fastq -t {threads} -o {params.qcdir} -a {input.adapters_tsv} {input.fq1} {input.fq2}
        
#         seqtk mergepe {input.fq1} {input.fq2} | \
#         bbduk.sh -Xmx4g -t={threads} in=stdin.fq out=stdout.fq interleaved=t ref={input.adapters_fa} \
#         ftm=5 minlen=25 qtrim=rl trimq=10 ktrim=r k=23 mink=11 hdist=1 tpe tbo | \
#         bbduk.sh -Xmx4g -t={threads} in=stdin.fq out={output.fq12} interleaved=t \
#         outm1={params.synt_fq1} outm2={params.synt_fq2} ref={input.refsynt_fa} k=31 hdist=1 stats={params.stats_synt}

#         fastqc --extract --delete -t {threads} -o {params.qcdir} -a {input.adapters_tsv} {output.fq12}

#         zip {output.qc_zip} {params.qcdir}/*
#         rm -r {params.qcdir}
#         """
