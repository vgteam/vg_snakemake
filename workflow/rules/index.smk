rule convert_gfa_to_gbz:
    input: config['gfa']
    output: 'results/pg/{graph}.gbz'
    threads: 8
    benchmark: 'benchmark/results/pg/{graph}.convert_gfa_to_gbz.benchmark.tsv'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    shell: "vg gbwt --num-jobs {threads} --gbz-format -g {output} -G {input}"

rule index_distance:
    input: '{graph}.gbz'
    output: '{graph}.dist'
    benchmark: 'benchmark/{graph}.index_distance.benchmark.tsv'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    shell: "vg index -j {output} {input}"

rule index_minimizer:
    input:
        gbz='{graph}.gbz',
        dist='{graph}.dist'
    threads: 8
    output: '{graph}.min'
    benchmark: 'benchmark/{graph}.index_minimizer.benchmark.tsv'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    shell: "vg minimizer -t {threads} -d {input.dist} -o {output} {input.gbz}"

rule index_r:
    input: "{graph}.gbz"
    output: "{graph}.ri"
    threads: 8
    benchmark: 'benchmark/{graph}.index_r.benchmark.tsv'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    shell: "vg gbwt -p --num-threads {threads} -r {output} -Z {input}"

rule index_haplotype_kmers:
    input:
        gbz='{graph}.gbz',
        dist='{graph}.dist',
        r='{graph}.ri'
    output: "{graph}.hapl"
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
    input: "{graph}.gbz"
    output: "{graph}.snarls"
    threads: 2
    benchmark: 'benchmark/{graph}.index_snarls.benchmark.tsv'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    shell: "vg snarls -T -t {threads} {input} > {output}"

rule extract_ref_fasta:
    input:
        gbz="{graph}.gbz",
        paths_list=config['ref_paths_list']
    output: "{graph}.ref.fa"
    benchmark: 'benchmark/{graph}.extract_ref_fasta.benchmark.tsv'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    shell: "vg paths --extract-fasta -p {input.paths_list} --xg {input.gbz} > {output}"

rule index_fasta:
    input: "{file}.fa"
    output: "{file}.fa.fai"
    benchmark: 'benchmark/{file}.index_fasta.benchmark.tsv'
    container: "docker://quay.io/biocontainers/samtools:1.18--hd87286a_0"
    shell: "samtools faidx {input}"

rule index_bam:
    input: "{file}.bam"
    output: "{file}.bam.bai"
    benchmark: 'benchmark/{file}.index_bam.benchmark.tsv'
    container: "docker://quay.io/biocontainers/samtools:1.18--hd87286a_0"
    shell: "samtools index {input} {output}"
