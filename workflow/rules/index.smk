rule convert_gfa_to_gbz:
    input: config['gfa']
    output: 'results/pg/{graph}.gbz'
    threads: 8
    benchmark: 'benchmark/{graph}.convert_gfa_to_gbz.benchmark.tsv'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    shell: "vg gbwt --num-jobs {threads} --gbz-format -g {output} -G {input}"

rule index_r_fullpg:
    input: 'results/pg/{graph}.gbz'
    output: 'results/pg/{graph}.ri'
    threads: 8
    benchmark: 'benchmark/{graph}.index_r.benchmark.tsv'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    shell: "vg gbwt -p --num-threads {threads} -r {output} -Z {input}"

rule index_distance_fullpg:
    input: 'results/pg/{graph}.gbz'
    output: 'results/pg/{graph}.dist'
    benchmark: 'benchmark/{graph}.index_distance.benchmark.tsv'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    shell: "vg index -j {output} {input}"

rule index_haplotype_kmers:
    input:
        gbz='results/pg/{graph}.gbz',
        dist='results/pg/{graph}.dist',
        r='results/pg/{graph}.ri'
    output: "results/pg/{graph}.hapl"
    threads: 8
    benchmark: 'benchmark/{graph}.index_haplotype_kmers.benchmark.tsv'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    shell:
        """
        vg haplotypes -v 2 \
        --kmer-length 29 --window-length 11 \
        --subchain-length 10000 -t {threads} \
        -d {input.dist} -r {input.r} -H {output} {input.gbz}
        """

rule index_snarls:
    input: 'results/pg/{graph}.gbz'
    output: 'results/pg/{graph}.snarls'
    threads: 2
    benchmark: 'benchmark/{graph}.index_snarls.benchmark.tsv'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    shell: "vg snarls -T -t {threads} {input} > {output}"

rule index_distance:
    input: "results/{sample}/{graph}.sample_pg.{sample}.gbz"
    output: temp("results/{sample}/{graph}.sample_pg.{sample}.dist")
    benchmark: 'benchmark/{sample}.{graph}.index_distance.benchmark.tsv'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    shell: "vg index -j {output} {input}"

rule index_minimizer:
    input:
        gbz="results/{sample}/{graph}.sample_pg.{sample}.gbz",
        dist="results/{sample}/{graph}.sample_pg.{sample}.dist"
    threads: 8
    output: temp("results/{sample}/{graph}.sample_pg.{sample}.min")
    benchmark: 'benchmark/{sample}.{graph}.index_minimizer.benchmark.tsv'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    shell: "vg minimizer -t {threads} -d {input.dist} -o {output} {input.gbz}"

rule extract_ref_fasta:
    input:
        gbz=getgbz(),
        paths_list=config['ref_paths_list']
    output: "results/pg/{graph}.ref.fa"
    benchmark: 'benchmark/{graph}.extract_ref_fasta.benchmark.tsv'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    params:
        seqn_prefix=config['seqn_prefix'],
        temp_paths='temp.path_name.{graph}.txt'
    shell:
        """
        rm -f {output} {params.temp_paths}
        for CHR in `cat {input.paths_list}`
        do
        echo $CHR > {params.temp_paths}
        vg paths --extract-fasta -p {params.temp_paths} --xg {input.gbz} | sed -e "s/>{params.seqn_prefix}/>/g" >> {output}
        done
        rm -f {params.temp_paths}
        """

rule index_fasta:
    input: getref()
    output:
        ref_idx="results/pg/{graph}.ref.fa.fai",
        dict="results/pg/{graph}.ref.dict"
    container: "docker://quay.io/cmarkello/samtools_picard@sha256:e484603c61e1753c349410f0901a7ba43a2e5eb1c6ce9a240b7f737bba661eb4"
    shell:
        """
        samtools faidx {input}
        # Save a reference copy by making the dict now
        java -jar /usr/picard/picard.jar CreateSequenceDictionary \
          R={input} \
          O={output.dict}
        """

rule index_bam:
    input: "results/{sample}/{sample}.{graph}.bam"
    output: temp("results/{sample}/{sample}.{graph}.bam.bai")
    benchmark: 'benchmark/{sample}.{graph}.index_bam.benchmark.tsv'
    container: "docker://quay.io/biocontainers/samtools:1.18--hd87286a_0"
    shell: "samtools index {input} {output}"

rule index_bam_tmp:
    input: "results/{sample}/temp_{sample}.{graph}.bam"
    output: temp("results/{sample}/temp_{sample}.{graph}.bam.bai")
    benchmark: 'benchmark/{sample}.{graph}.index_bam.benchmark.tsv'
    container: "docker://quay.io/biocontainers/samtools:1.18--hd87286a_0"
    shell: "samtools index {input} {output}"
