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
    output: "results/{sample}/{graph}.sample_pg.{sample}.dist"
    benchmark: 'benchmark/{sample}.{graph}.index_distance.benchmark.tsv'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    shell: "vg index -j {output} {input}"

rule index_minimizer:
    input:
        gbz="results/{sample}/{graph}.sample_pg.{sample}.gbz",
        dist="results/{sample}/{graph}.sample_pg.{sample}.dist"
    threads: 8
    output: "results/{sample}/{graph}.sample_pg.{sample}.min"
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
    shell: "vg paths --extract-fasta -p {input.paths_list} --xg {input.gbz} > {output}"

rule index_fasta:
    input: "results/pg/{graph}.ref.fa"
    output: "results/pg/{graph}.ref.fa.fai"
    benchmark: 'benchmark/{graph}.index_fasta.benchmark.tsv'
    container: "docker://quay.io/biocontainers/samtools:1.18--hd87286a_0"
    shell: "samtools faidx {input}"

rule index_bam:
    input: "results/{sample}/{sample}.{graph}.bam"
    output: "results/{sample}/{sample}.{graph}.bam.bai"
    benchmark: 'benchmark/{sample}.{graph}.index_bam.benchmark.tsv'
    container: "docker://quay.io/biocontainers/samtools:1.18--hd87286a_0"
    shell: "samtools index {input} {output}"
