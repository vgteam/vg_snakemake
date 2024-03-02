if len(config['refsynt_fa']) > 0 and len(config['adapters_fa']) > 0 and len(config['adapters_tsv']) > 0:
    rule count_kmer_in_reads:
        input: "results/{sample}/{sample}.trimmed.fq.gz"
        output: temp("results/{sample}/{sample}.kff")
        params:
            wdir="temp_kmc_{sample}",
            filel="temp_kmc_filelist_{sample}.txt"
        threads: 8
        benchmark: 'benchmark/{sample}.count_kmer_in_reads.benchmark.tsv'
        container: "docker://quay.io/biocontainers/kmc:3.2.1--hf1761c0_2"
        shell:
            """
            echo {input} > {params.filel}
            rm -rf {params.wdir}
            mkdir -p {params.wdir}
            kmc -k29 -m64 -okff -t{threads} @{params.filel} {params.wdir}/out {params.wdir}
            mv {params.wdir}/out.kff {output}
            rm -r {params.filel} {params.wdir}
            """

    # if reads were trimmed, align from the trimmed interleaved fastq
    rule map_short_reads_giraffe:
        input: 
            fq12="results/{sample}/{sample}.trimmed.fq.gz",
            gbz="results/{sample}/{graph}.sample_pg.{sample}.gbz",
            dist='results/{sample}/{graph}.sample_pg.{sample}.dist',
            min='results/{sample}/{graph}.sample_pg.{sample}.min'
        output: temp("results/{sample}/{sample}.{graph}.gaf.gz")
        threads: 8
        priority: 2
        benchmark: 'benchmark/{sample}.{graph}.map_short_reads_giraffe.benchmark.tsv'
        container: "docker://quay.io/vgteam/vg:v1.52.0"
        shell:
            """
            vg giraffe --progress \
            --sample "{wildcards.sample}" \
            --output-format gaf \
            -f {input.fq12} -i \
            -Z {input.gbz} \
            -d {input.dist} \
            -m {input.min} \
            -t {threads} | gzip > {output}
            """
else:
    rule count_kmer_in_reads:
        input:
            fq1=getfq1,
            fq2=getfq2
        output: temp("results/{sample}/{sample}.kff")
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
            rm -rf {params.wdir}
            mkdir -p {params.wdir}
            kmc -k29 -m64 -okff -t{threads} @{params.filel} {params.wdir}/out {params.wdir}
            mv {params.wdir}/out.kff {output}
            rm -r {params.filel} {params.wdir}
            """

    # otherwise, align from the input fastq pair
    rule map_short_reads_giraffe:
        input: 
            fq1=getfq1,
            fq2=getfq2,
            gbz="results/{sample}/{graph}.sample_pg.{sample}.gbz",
            dist='results/{sample}/{graph}.sample_pg.{sample}.dist',
            min='results/{sample}/{graph}.sample_pg.{sample}.min'
        output: temp("results/{sample}/{sample}.{graph}.gaf.gz")
        threads: 8
        priority: 2
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
    output: temp("results/{sample}/{graph}.sample_pg.{sample}.gbz")
    threads: 8
    priority: 1
    benchmark: 'benchmark/{sample}.{graph}.sample_haplotypes.benchmark.tsv'
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
        gaf="results/{sample}/{sample}.{graph}.gaf.gz",
        ref=getref(),
        ref_idx=getrefidx(),
        rename_script="workflow/scripts/rename_bam_stream.py"
    output: temp("results/{sample}/{sample}.{graph}.surj.bam")
    priority: 3
    params:
        surj_threads=lambda wildcards, threads: max(1, int(threads/2)) if threads < 8 else threads - 4,
        sort_threads=lambda wildcards, threads: max(1, int(threads/2)) if threads < 8 else 4,
        sort_dir="temp_bam_sort_{sample}_{graph}",
        seqn_prefix=config['seqn_prefix']
    threads: 8
    benchmark: 'benchmark/{sample}.{graph}.surject_reads.benchmark.tsv'
    container: "docker://quay.io/jmonlong/freebayes-samtools-vg:1.2.0_1.10.0_1.53.0"
    shell:
        """
        rm -rf {params.sort_dir}
        mkdir -p {params.sort_dir}

        vg surject \
        -F {input.paths_list} \
        -x {input.gbz} \
        -t {params.surj_threads} \
        --sam-output --gaf-input \
        --sample {wildcards.sample} \
        --read-group "ID:1 LB:lib1 SM:{wildcards.sample} PL:illumina PU:unit1" \
        --prune-low-cplx --interleaved --max-frag-len 3000 \
        {input.gaf} | \
        python {input.rename_script} -f {input.ref_idx} -p "{params.seqn_prefix}" | \
        bamleftalign --fasta-reference {input.ref} --compressed | \
        samtools sort -m 4G --threads {params.sort_threads} -T {params.sort_dir}/temp \
        -O BAM > {output}

        rm -rf {params.sort_dir}
        """

use rule surject_reads as surject_reads_tmp with :
    output: temp("results/{sample}/temp_{sample}.{graph}.surj.bam")

rule prepare_target_regions:
    input: 
        ref=getref(),
        ref_idx=getrefidx(),
        ref_dict=getrefdict(),
        bam="results/{sample}/temp_{sample}.{graph}.surj.bam",
        bai="results/{sample}/temp_{sample}.{graph}.surj.bam.bai"
    output: temp("results/{sample}/{sample}.{graph}.realn_targets.bed")
    container: "docker://quay.io/jmonlong/gatk-bedtools:3.8.1_2.21.0"
    threads: 8
    priority: 4
    benchmark: "benchmark/{sample}.{graph}.prepare_target_regions.benchmark.tsv"
    params:
        inter="{sample}.{graph}.forIndelRealigner.intervals",
        bed="{sample}.{graph}.forIndelRealigner.bed"
    shell:
        """
        java -Xmx16G -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \
          --remove_program_records \
          -drf DuplicateRead \
          --disable_bam_indexing \
          -nt {threads} \
          -R {input.ref} \
          -I {input.bam} \
          --out {params.inter}

        awk -F '[:-]' 'BEGIN {{ OFS = "\t" }} {{ if( $3 == "") {{ print $1, $2-1, $2 }} else {{ print $1, $2-1, $3}}}}' {params.inter} > {params.bed}

        bedtools slop -i {params.bed} -g "{input.ref_idx}" -b "160" > {output}
        rm -f {params.bed} {params.inter}
        """

rule realign_bam:
    input:
        ref=getref(),
        ref_idx=getrefidx(),
        target_bed="results/{sample}/{sample}.{graph}.realn_targets.bed",
        bam="results/{sample}/temp_{sample}.{graph}.surj.bam",
        bai="results/{sample}/temp_{sample}.{graph}.surj.bam.bai"
    output: temp("results/{sample}/{sample}.{graph}.surj_realn.bam")
    params:
        tmpdir="temp.realign_bam.{sample}.{graph}"
    benchmark: "benchmark/{sample}.{graph}.realign_bam.benchmark.tsv"
    log: "logs/{sample}.{graph}.realign_bam.log"
    priority: 5
    container: 'docker://quay.io/adamnovak/dceoy-abra2@sha256:43d09d1c10220cfeab09e2763c2c5257884fa4457bcaa224f4e3796a28a24bba'
    threads: 8
    shell:
        """
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir}
        java -Xmx28G -jar /opt/abra2/abra2.jar \
          --targets {input.target_bed} \
          --in {input.bam} \
          --out {output} \
          --tmpdir {params.tmpdir} \
          --ref {input.ref} \
          --threads {threads} 2> {log}
        rm -rf {params.tmpdir}
        """

if len(config['refsynt_fa']) > 0 and len(config['adapters_fa']) > 0 and len(config['adapters_tsv']) > 0:
    rule extract_unmapped_reads:
        input:
            fq="results/{sample}/{sample}.trimmed.fq.gz",
            gaf="results/{sample}/{sample}.{graph}.gaf.gz"
        output: temp("results/{sample}/{sample}.{graph}.unmapped.fq.gz")
        params:
            reads="temp.{sample}.{graph}.unmapped.reads.txt"
        threads: 1
        benchmark: 'benchmark/{sample}.{graph}.extract_unmapped_reads.benchmark.tsv'    
        shell:
            """
            zcat {input.gaf} | awk '{{if($3=="*"){{print $1}}}}' | uniq > {params.reads}
            seqtk subseq {input.fq} {params.reads} | gzip > {output}
            rm {params.reads}
            """
else:
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
    

rule trim_fastq:
    input:
        fq1=getfq1,
        fq2=getfq2,
        refsynt_fa=config['refsynt_fa'],
        adpt_fa=config['adapters_fa']
    output:
        fq12=temp("results/{sample}/{sample}.trimmed.fq.gz"),
        stats_synt=temp("results/{sample}/{sample}.statsSYNT.txt"),
        synt_fq1=temp("results/{sample}/{sample}.synt.1.fq.gz"),
        synt_fq2=temp("results/{sample}/{sample}.synt.2.fq.gz")
    benchmark: 'benchmark/{sample}.trim_fastq.benchmark.tsv'    
    threads: 4
    shell:
        """
        seqtk mergepe {input.fq1} {input.fq2} | \
        bbduk.sh -Xmx7g -t={threads} in=stdin.fq out=stdout.fq interleaved=t ref={input.adpt_fa} \
        ftm=5 minlen=25 qtrim=rl trimq=10 ktrim=r k=23 mink=11 hdist=1 tpe tbo | \
        bbduk.sh -Xmx7g -t={threads} in=stdin.fq out={output.fq12} interleaved=t \
        outm1={output.synt_fq1} outm2={output.synt_fq2} ref={input.refsynt_fa} k=31 hdist=1 stats={output.stats_synt}
        """

rule qc_fastq:
    input:
        fq1=getfq1,
        fq2=getfq2,
        fqt="results/{sample}/{sample}.trimmed.fq.gz",
        adpt_tsv=config['adapters_tsv'],
        stats_synt="results/{sample}/{sample}.statsSYNT.txt",
        synt_fq1="results/{sample}/{sample}.synt.1.fq.gz",
        synt_fq2="results/{sample}/{sample}.synt.2.fq.gz"
    output: temp("results/{sample}/fastq_qc.{sample}.zip")
    params:
        qcdir="temp_fastqc.{sample}",
        stats_synt="temp_fastqc.{sample}/{sample}.statsSYNT.txt",
        synt_fq1="temp_fastqc.{sample}/{sample}.synt.1.fq.gz",
        synt_fq2="temp_fastqc.{sample}/{sample}.synt.2.fq.gz"
    benchmark: 'benchmark/{sample}.qc_fastq.benchmark.tsv'    
    threads: 4
    shell:
        """
        rm -rf {params.qcdir}
        mkdir -p {params.qcdir}
        fastqc --extract --delete -f fastq -t {threads} -o {params.qcdir} -a {input.adpt_tsv} {input.fq1} {input.fq2} {input.fqt}
        cp {input.stats_synt} {input.synt_fq1} {input.synt_fq2} {params.qcdir}/
        zip {output.qc_zip} {params.qcdir}/*
        rm -r {params.qcdir}
        """
