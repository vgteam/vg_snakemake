# Run with, for example:
# > snakemake --configfile config.yaml --resources mem_mb=200000 --cores 16 construct_all
# Config can also be overwritten in the command line, for example:
# > snakemake --configfile config.yaml --config samples="HG00514" --resources mem_mb=100000 --cores 16 genotype

##
## Config
##

# chromosome names
CHRS=[config['chr_prefix'] + ii for ii in config['chrs'].split()]

# config values used a lot directly within shell commands
SROOT=config['s3root']
GRAPH=config['graph']
# parse config values
SAMPLES=config['samples'].split()
MAPPER=config['mapper']
if MAPPER == 'giraffe':
    MAPPER = 'giraffe{}k{}w{}N'.format(config['mink'], config['minw'], config['covern'])

##
## Main rules
##

# indexes used by the default mapper and variant caller
rule construct_all:
    input: expand('{graph}.{ext}', graph=GRAPH, ext=['xg', 'snarls', 'gcsa'])

# indexes used by the giraffe mapper and variant caller
rule construct_all_giraffe:
    input:
        expand('{graph}.{ext}', graph=GRAPH, ext=['xg', 'snarls', 'dist']),
        expand('{graph}.k{k}.w{w}.N{n}.min', graph=GRAPH, k=config['mink'],
               w=config['minw'], n=config['covern'])

# genotype variants in the VCF used to build the graph
rule genotype:
    input: expand('{sample}/{sample}-{graph}.{map}.q5.gt.vcf.gz', sample=SAMPLES, graph=GRAPH, map=MAPPER)

# call variants in the graph
rule call:
    input: expand('{sample}/{sample}-{graph}.{map}.q5.call.vcf.gz', sample=SAMPLES, graph=GRAPH, map=MAPPER)

# augment and call variants in the graph
rule augcall:
    input: expand('{sample}/{sample}-{graph}.{map}.aug.q5.call.vcf.gz', sample=SAMPLES, graph=GRAPH, map=MAPPER)

# remove ALL files (only use to start again from scratch)
rule cleanall:
    params:
        ff=expand('{graph}.{ext}', graph=GRAPH, ext=['xg', 'snarls', 'gcsa', 'gcsa.lcp', 'ids.mapping', 'dist', 'trivial.snarls']) +
        expand('{graph}-{chrom}{prune}.vg', graph=GRAPH, chrom=CHRS, prune=['', '-pruned']) +
        expand('{graph}.k{k}.w{w}.N{n}.min', graph=GRAPH, k=config['mink'],
               w=config['minw'], n=config['covern']) +
        expand('{graph}.N{n}.{ext}', graph=GRAPH, n=config['covern'], ext=['gbwt', 'gg']) +
        expand('{sample}/{sample}-{graph}.{map}.*', sample=SAMPLES, graph=GRAPH, map=MAPPER) +
        expand('{sample}/read_chunks', sample=SAMPLES, graph=GRAPH, map=MAPPER)
    shell:
        "rm -rf {params.ff}"

#
# Graph construction
#

# make a .vg for a chromosome
rule construct_chr:
    input:
        ref='{genome}.fa',
        vcf='{svs}.vcf.gz',
        vcfidx='{svs}.vcf.gz.tbi'
    output: '{genome}-{svs}-{chr}.vg'
    benchmark: 'benchmarks/{genome}-{svs}-{chr}.construct.benchmark.txt'
    log: "logs/{genome}-{svs}-{chr}.construct.log.txt"
    threads: 1
    run:
        shell('(vg construct -t {threads} -r {input.ref} -v {input.vcf} -a -f -S -R {wildcards.chr} -C | vg ids --sort - > {output}) 2> {log}')
        if config['s3save']:
            shell("aws s3 cp --quiet {output} {SROOT}/")

# join ids across the multiple .vg files.
rule join_ids:
    input: expand('{{genome}}-{{svs}}-{chr}.vg', chr=CHRS)
    output:
        mapping='{genome}-{svs}.ids.mapping'
    threads: 1
    benchmark: 'benchmarks/{genome}-{svs}-joinids.benchmark.txt'
    log: 'logs/{genome}-{svs}-joinids.log.txt'
    run:
        shell('vg ids --join --mapping {output.mapping} {input} 2> {log}')
        if config['s3save']:
            shell("for ff in {input}; do aws s3 cp --quiet $ff {SROOT}/; done")
            shell("aws s3 cp --quiet {output.mapping} {SROOT}/")

# make xg index containing the alts paths
rule index_xg:
    input:
        vg=expand('{{genome}}-{{svs}}-{chr}.vg', chr=CHRS),
        mapping='{genome}-{svs}.ids.mapping'   
    output: '{genome}-{svs}.xg'
    threads: config['cores_xg']
    resources:
        mem_mb=config['mem_xg']
    benchmark: 'benchmarks/{genome}-{svs}-xg.benchmark.txt'
    log: 'logs/{genome}-{svs}-xg.log.txt'
    run:
        shell('vg index -L -t {threads} -x {output} {input.vg} 2> {log}')
        if config['s3save']:
            shell("aws s3 cp --quiet {output} {SROOT}/")

# prepare the snarls index
rule index_snarls:
    input: '{genome}-{svs}.xg'
    output: '{genome}-{svs}.snarls'
    threads: config['cores_snarls']
    resources:
        mem_mb=config['mem_snarls']
    benchmark: 'benchmarks/{genome}-{svs}-snarls.benchmark.txt'
    log: 'logs/{genome}-{svs}-snarls.log.txt'
    run:
        shell('vg snarls -t {threads} {input} > {output} 2> {log}')
        if config['s3save']:
            shell("aws s3 cp --quiet {output} {SROOT}/")

# pruning is used before constructing the GSCA index
rule prune_vg:
    input:
        vg='{genome}-{svs}-{chr}.vg',
        mapping='{genome}-{svs}.ids.mapping'
    output: '{genome}-{svs}-{chr}-pruned.vg'
    threads: 1
    benchmark: 'benchmarks/{genome}-{svs}-{chr}.prune.benchmark.txt'
    log: 'logs/{genome}-{svs}-{chr}.prune.log.txt'
    run:
        shell('vg prune -t {threads} -M 32 --restore-paths {input.vg} > {output} 2> {log}')
        if config['s3save']:
            shell("aws s3 cp --quiet {output} {SROOT}/")        

# GCSA index construction from the pruned graphs
rule index_gcsa:
    input: expand('{{genome}}-{{svs}}-{chr}-pruned.vg', chr=CHRS)
    output:
        gcsa='{genome}-{svs}.gcsa',
        gcsalcp='{genome}-{svs}.gcsa.lcp'
    threads: config['cores_gcsa']
    resources:
        mem_mb=config['mem_gcsa']
    benchmark: 'benchmarks/{genome}-{svs}-gcsa.benchmark.txt'
    log: 'logs/{genome}-{svs}-gcsa.log.txt'
    params:
        tmp_dir="temp_gsca_{genome}_{svs}"
    run:
        shell('mkdir -p {params.tmp_dir}')
        shell('vg index --temp-dir {params.tmp_dir} -p -t {threads} -g {output.gcsa} {input} 2> {log}')
        shell('rm -r {params.tmp_dir}')
        if config['s3save']:
            shell("aws s3 cp --quiet {output.gcsa} {SROOT}/")
            shell("aws s3 cp --quiet {output.gcsalcp} {SROOT}/")

## GBWT with greedy path cover
rule index_gbwt_greedy:
    input: '{graph}.xg'
    output:
        gg='{graph}.N{n}.gg',
        gbwt='{graph}.N{n}.gbwt'
    threads: 1
    resources:
        mem_mb=config['mem_gbwt']
    benchmark: 'benchmarks/{graph}-gbwt-N{n}.benchmark.txt'
    log: 'logs/{graph}-gbwt-N{n}.log.txt'
    run:
        shell('vg gbwt -n {wildcards.n} -g {output.gg} -o {output.gbwt} -x {input} -P 2> {log}')
        if config['s3save']:
            shell("aws s3 cp --quiet {output.gbwt} {SROOT}/")
            shell("aws s3 cp --quiet {output.gg} {SROOT}/")

rule index_minimizer:
    input:
        xg='{graph}.xg',
        gbwt='{graph}.N{n}.gbwt'
    output: '{graph}.k{k}.w{w}.N{n}.min'
    threads: config['cores_minimizer_index']
    resources:
        mem_mb=config['mem_minimizer_index']
    benchmark: 'benchmarks/{graph}-minimizer-k{k}-w{w}-N{n}.benchmark.txt'
    log: 'logs/{graph}-minimizer-k{k}-w{w}-N{n}.log.txt'
    run:
        shell('vg minimizer -k {wildcards.k} -w {wildcards.w} -t {threads} -i {output} -g {input.gbwt} {input.xg} 2> {log}')
        if config['s3save']:
            shell("aws s3 cp --quiet {output} {SROOT}/")

rule index_trivial_snarls:
    input: '{genome}-{svs}.xg'
    output: '{genome}-{svs}.trivial.snarls'
    threads: config['cores_snarls']
    resources:
        mem_mb=config['mem_snarls']
    benchmark: 'benchmarks/{genome}-{svs}-trivialsnarls.benchmark.txt'
    log: 'logs/{genome}-{svs}-trivialsnarls.log.txt'
    run:
        shell('vg snarls -t {threads} --include-trivial {input} > {output} 2> {log}')
        if config['s3save']:
            shell("aws s3 cp --quiet {output} {SROOT}/")

rule index_distance:
    input:
        xg='{graph}.xg',
        snarls='{graph}.trivial.snarls'
    output: '{graph}.dist'
    threads: config['cores_dist_index']
    resources:
        mem_mb=config['mem_dist_index']
    benchmark: 'benchmarks/{graph}-distance.benchmark.txt'
    log: 'logs/{graph}-distance.log.txt'
    run:
        shell('vg index -t {threads} -j {output} -x {input.xg} -s {input.snarls} 2> {log}')
        if config['s3save']:
            shell("aws s3 cp --quiet {output} {SROOT}/")

#
# Map reads from a sample and call variants
#

# Eventually split the reads into chunk to map in parallel
if config['nb_split_reads'] > 0:
    # split fastq into chunks
    checkpoint split_reads_1:
        input: '{sample}/{sample}_1.fastq.gz'
        output:
            dir=directory('{sample}/read_chunks/1')
        threads: config['cores_split_reads']
        benchmark: 'benchmarks/{sample}-1-splitreads.benchmark.txt'
        run:
            NLINES=config['nb_split_reads'] * 4
            shell('mkdir -p {output.dir}')
            shell("gzip -cd {input} | split -d -l {NLINES} --filter='pigz -p {threads} > ${{FILE}}.fastq.gz' - \"{output.dir}/{wildcards.sample}_1.part\"")
    checkpoint split_reads_2:
        input: '{sample}/{sample}_2.fastq.gz'
        output:
            dir=directory('{sample}/read_chunks/2')
        threads: config['cores_split_reads']
        benchmark: 'benchmarks/{sample}-2-splitreads.benchmark.txt'
        run:
            NLINES=config['nb_split_reads'] * 4
            shell('mkdir -p {output.dir}')
            shell("gzip -cd {input} | split -d -l {NLINES} --filter='pigz -p {threads} > ${{FILE}}.fastq.gz' - \"{output.dir}/{wildcards.sample}_2.part\"")
    # set the input paths for the mapping rules
    read1_in = '{sample}/read_chunks/1/{sample}_1.part{part}.fastq.gz'
    read2_in = '{sample}/read_chunks/2/{sample}_2.part{part}.fastq.gz'
    map_lab = '{sample}.part{part}'
    map_out = '{sample}/read_chunks/{sample}-{genome}-{svs}.part{part}'
else:
    # set the input paths for the mapping rules (no chunks)
    read1_in = '{sample}/{sample}_1.fastq.gz'
    read2_in = '{sample}/{sample}_2.fastq.gz'
    map_lab = '{sample}'
    map_out = '{sample}/{sample}-{genome}-{svs}'
    
# map reads to the graph
rule map:
    input:
        r1=read1_in,
        r2=read2_in,
        xg='{genome}-{svs}.xg',
        gcsa='{genome}-{svs}.gcsa',
        gcsalcp='{genome}-{svs}.gcsa.lcp'
    output: map_out + '.map.gam'
    threads: config['cores_map']
    resources:
        mem_mb=config['mem_map']
    benchmark: 'benchmarks/' + map_lab + '.{genome}.{svs}.map.benchmark.txt'
    log: 'logs/' + map_lab + '-{genome}-{svs}-map.log.txt'
    run:
        shell("vg map -t {threads} -x {input.xg} -g {input.gcsa} -f {input.r1} -f {input.r2} > {output} 2> {log}")
        if config['s3save']:
            shell("aws s3 cp --quiet {output} {SROOT}/{output}")

# map reads to the graph using mpmap in single-path mode
rule map_mpmap:
    input:
        r1=read1_in,
        r2=read2_in,
        xg='{genome}-{svs}.xg',
        gcsa='{genome}-{svs}.gcsa',
        gcsalcp='{genome}-{svs}.gcsa.lcp'
    output: map_out + '.mpmap.gam'
    threads: config['cores_map']
    resources:
        mem_mb=config['mem_map']
    benchmark: 'benchmarks/' + map_lab + '.{genome}.{svs}.mpmap.benchmark.txt'
    log: 'logs/' + map_lab + '-{genome}-{svs}-mpmap.log.txt'
    run:
        shell("vg mpmap -S -t {threads} -x {input.xg} -g {input.gcsa} -f {input.r1} -f {input.r2} > {output} 2> {log}")
        if config['s3save']:
            shell("aws s3 cp --quiet {output} {SROOT}/{output}")

# map reads to the graph using giraffe
rule map_giraffe:
    input:
        r1=read1_in,
        r2=read2_in,
        xg='{genome}-{svs}.xg',
        min='{genome}-{svs}.k{k}.w{w}.N{n}.min',
        dist='{genome}-{svs}.dist',
        gbwt='{genome}-{svs}.N{n}.gbwt'
    output: map_out + '.giraffe{k}k{w}w{n}N.gam'
    threads: config['cores_map']
    resources:
        mem_mb=config['mem_map']
    benchmark: 'benchmarks/' + map_lab + '.{genome}.{svs}.giraffe{k}k{w}w{n}N.benchmark.txt'
    log: 'logs/' + map_lab + '-{genome}-{svs}-giraffe{k}k{w}w{n}N.log.txt'
    run:
        shell("vg giraffe -p -t {threads} -m {input.min} -d {input.dist} --gbwt-name {input.gbwt} -x {input.xg} -N {wildcards.sample} -f {input.r1} -f {input.r2} > {output} 2> {log}")
        if config['s3save']:
            shell("aws s3 cp --quiet {output} {SROOT}/{output}")


# If the reads were split, they need to be merged back after alignment
if config['nb_split_reads'] > 0:
    # merge aligned reads
    def aggregate_reads(wildcards):
        checkpoint_output = checkpoints.split_reads_1.get(**wildcards).output[0]
        checkpoint_output2 = checkpoints.split_reads_2.get(**wildcards).output[0]
        return expand("{sample}/read_chunks/{sample}-{graph}.part{part}.{map}.gam",
                      sample=wildcards.sample, graph=wildcards.graph, map=wildcards.map,
                      part=glob_wildcards(os.path.join(checkpoint_output, wildcards.sample + "_1.part{part}.fastq.gz")).part)
    rule merge_gam:
        input: aggregate_reads
        output: '{sample}/{sample}-{graph}.{map}.gam'
        benchmark: 'benchmarks/{sample}-{graph}-{map}-mergegam.benchmark.txt'
        run:
            shell("cat {input} > {output}")
            if config['s3save']:
                shell("aws s3 cp --quiet {output} {SROOT}/{wildcards.sample}/")


# compute packed coverage from aligned reads
rule pack:
    input:
        gam='{sample}/{sample}-{graph}.{map}.gam',
        xg='{graph}.xg'
    output: '{sample}/{sample}-{graph}.{map}.q{minq}.pack'
    threads: config['cores_pack']
    resources:
        mem_mb=config['mem_pack']
    benchmark: 'benchmarks/{sample}-{graph}-{map}-q{minq}-pack.benchmark.txt'
    log: 'logs/{sample}-{graph}-{map}-q{minq}-pack.log.txt'
    run:
        shell("vg pack -x {input.xg} -g {input.gam} -Q {wildcards.minq} -t {threads} -o {output} 2> {log}")
        if config['s3save']:
            shell("aws s3 cp --quiet {output} {SROOT}/{wildcards.sample}/")

# call variants from the packed read coverage
rule call_novcf:
    input:
        pack='{sample}/{sample}-{graph}.{map}.q{minq}.pack',
        xg='{graph}.xg',
        snarls='{graph}.snarls'
    output:
        vcf='{sample}/{sample}-{graph}.{map}.q{minq}.call.vcf.gz',
        idx='{sample}/{sample}-{graph}.{map}.q{minq}.call.vcf.gz.tbi'
    params:
        tmp_raw_vcf="{sample}-{graph}-q{minq}_calltemp_raw.vcf"
    threads: config['cores_call']
    resources:
        mem_mb=config['mem_call']
    benchmark: 'benchmarks/{sample}-{graph}-{map}-q{minq}-call.benchmark.txt'
    log: 'logs/{sample}-{graph}-{map}-q{minq}-call.log.txt'
    run:
        shell("vg call -k {input.pack} -t {threads} -s {wildcards.sample} --snarls {input.snarls} {input.xg} > {params.tmp_raw_vcf} 2> {log}")
        shell("bcftools sort {params.tmp_raw_vcf} | bgzip > {output.vcf}")
        shell("tabix -f -p vcf {output.vcf}")
        shell('rm {params.tmp_raw_vcf}')
        if config['s3save']:
            shell("aws s3 cp --quiet {output.vcf} {SROOT}/{wildcards.sample}/")
            shell("aws s3 cp --quiet {output.idx} {SROOT}/{wildcards.sample}/")

# genotype variants from the packed read coverage
rule call_vcf:
    input:
        pack='{sample}/{sample}-{genome}-{svs}.{map}.q{minq}.pack',
        xg='{genome}-{svs}.xg',
        snarls='{genome}-{svs}.snarls',
        vcf='{svs}.vcf.gz',
        vcftbi='{svs}.vcf.gz.tbi'
    output:
        vcf='{sample}/{sample}-{genome}-{svs}.{map}.q{minq}.gt.vcf.gz',
        idx='{sample}/{sample}-{genome}-{svs}.{map}.q{minq}.gt.vcf.gz.tbi'
    params:
        tmp_raw_vcf="{sample}-{genome}-{svs}-q{minq}_genotemp_raw.vcf"
    threads: config['cores_call']
    resources:
        mem_mb=config['mem_call']
    benchmark: 'benchmarks/{sample}-{genome}-{svs}-{map}-q{minq}-genotype.benchmark.txt'
    log: 'logs/{sample}-{genome}-{svs}-{map}-q{minq}-genotype.log.txt'
    run:
        shell("vg call -k {input.pack} -t {threads} -s {wildcards.sample} --snarls {input.snarls} -v {input.vcf} {input.xg} > {params.tmp_raw_vcf} 2> {log}")
        shell("bcftools view -e 'GT=\"0/0\" || GT=\"./.\"' {params.tmp_raw_vcf} | bcftools sort | bgzip > {output.vcf}")
        shell("tabix -f -p vcf {output.vcf}")
        shell('rm {params.tmp_raw_vcf}')
        if config['s3save']:
            shell("aws s3 cp --quiet {output.vcf} {SROOT}/{wildcards.sample}/")
            shell("aws s3 cp --quiet {output.idx} {SROOT}/{wildcards.sample}/")

#
# Augment graph and call variants
#

# convert a XG index to a PG index
rule pgconvert:
    input: '{graph}.xg'
    output: '{graph}.pg'
    threads: 1
    resources:
        mem_mb=config['mem_pgconvert']
    benchmark: 'benchmarks/{graph}-pgconvert.benchmark.txt'
    log: 'logs/{graph}-pgconvert.log.txt'
    run:
        shell('vg convert {input} -p > {output} 2> {log}')
        if config['s3save']:
            shell("aws s3 cp --quiet {output} {SROOT}/")

# augment a graph with aligned reads
rule augment:
    input:
        pg='{graph}.pg',
        gam='{sample}/{sample}-{graph}.{map}.gam'
    output:
        gam='{sample}/{sample}-{graph}.{map}.aug.gam',
        pg='{sample}/{sample}-{graph}.{map}.aug.pg'
    threads: config['cores_augment']
    resources:
        mem_mb=config['mem_augment']
    benchmark: 'benchmarks/{sample}-{graph}-{map}-augment.benchmark.txt'
    log: 'logs/{sample}-{graph}-{map}-augment.log.txt'             
    run:
        shell('vg augment {input.pg} {input.gam} -t {threads} -m 4 -q 5 -Q 5 -A {output.gam} > {output.pg} 2> {log}')
        if config['s3save']:
            shell("aws s3 cp --quiet {output.gam} {SROOT}/{wildcards.sample}/")
            shell("aws s3 cp --quiet {output.pg} {SROOT}/{wildcards.sample}/")

# prepare the snarls index for the augmented graph
rule index_snarls_aug:
    input: '{sample}/{sample}-{graph}.{map}.aug.pg'
    output: '{sample}/{sample}-{graph}.{map}.aug.snarls'
    threads: config['cores_snarls']
    resources:
        mem_mb=config['mem_snarls']
    benchmark: 'benchmarks/{sample}-{graph}-{map}-aug-snarls.benchmark.txt'
    log: 'logs/{sample}-{graph}-{map}-aug-snarls.log.txt'
    run:
        shell('vg snarls -t {threads} {input} > {output} 2> {log}')
        if config['s3save']:
            shell("aws s3 cp --quiet {output} {SROOT}/{wildcards.sample}/")

# compute packed coverage on augmented graph
rule pack_aug:
    input:
        gam='{sample}/{sample}-{graph}.{map}.gam',
        pg='{sample}/{sample}-{graph}.{map}.aug.pg'
    output: '{sample}/{sample}-{graph}.{map}.aug.q{minq}.pack'
    threads: config['cores_pack']
    resources:
        mem_mb=config['mem_pack']
    benchmark: 'benchmarks/{sample}-{graph}-{map}-aug-q{minq}-pack.benchmark.txt'
    log: 'logs/{sample}-{graph}-{map}-aug-q{minq}-pack.log.txt'
    run:
        shell("vg pack -x {input.pg} -g {input.gam} -Q {wildcards.minq} -t {threads} -o {output} 2> {log}")
        if config['s3save']:
            shell("aws s3 cp --quiet {output} {SROOT}/{wildcards.sample}/")

# call variants from the packed read coverage
rule call_aug:
    input:
        pack='{sample}/{sample}-{graph}.{map}.aug.q{minq}.pack',
        pg='{sample}/{sample}-{graph}.{map}.aug.pg',
        snarls='{sample}/{sample}-{graph}.{map}.aug.snarls'
    output:
        vcf='{sample}/{sample}-{graph}.{map}.aug.q{minq}.call.vcf.gz',
        idx='{sample}/{sample}-{graph}.{map}.aug.q{minq}.call.vcf.gz.tbi'
    params:
        tmp_raw_vcf="{sample}-{graph}-aug.q{minq}_calltemp_raw.vcf"
    threads: config['cores_call']
    resources:
        mem_mb=config['mem_call']
    benchmark: 'benchmarks/{sample}-{graph}-{map}-aug-q{minq}-call.benchmark.txt'
    log: 'logs/{sample}-{graph}-{map}-aug-q{minq}-call.log.txt'
    run:
        shell("vg call -k {input.pack} -t {threads} -s {wildcards.sample} --snarls {input.snarls} {input.xg} > {params.tmp_raw_vcf} 2> {log}")
        shell("bcftools sort {params.tmp_raw_vcf} | bgzip > {output.vcf}")
        shell("tabix -f -p vcf {output.vcf}")
        shell('rm {params.tmp_raw_vcf}')
        if config['s3save']:
            shell("aws s3 cp --quiet {output.vcf} {SROOT}/{wildcards.sample}/")
            shell("aws s3 cp --quiet {output.idx} {SROOT}/{wildcards.sample}/")
