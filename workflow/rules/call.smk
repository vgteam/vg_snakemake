rule pack:
    input:
        gbz=getgbz(),
        gaf="results/gaf/{sample}.{graph}.gaf.gz"
    output: 'results/vgcall/{sample}.{graph}.pack'
    threads: 8
    benchmark: 'benchmark/results/vgcall/{sample}.{graph}.pack.benchmark.tsv'
    log: 'logs/pack.{graph}.{sample}.log'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    shell: "vg pack -x {input.gbz} -a {input.gaf} -Q 5 -t {threads} -o {output} 2> {log}"

rule vgcall:
    input: 
        gbz=getgbz(),
        snarls=getsnarls(),
        paths_list=config['ref_paths_list'],
        pack='results/vgcall/{sample}.{graph}.pack'
    output: 'results/vcf/{sample}.{graph}.gt.vcf.gz'
    threads: 8
    benchmark: 'benchmark/results/vcf/{sample}.{graph}.vgcall.benchmark.tsv'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    shell:
        """
        RPATHS=""
        for RP in `cat {input.paths_list}`
        do
        RPATHS=`echo "-p $RP $RPATHS"`
        done
        
        vg call -t {threads} -r {input.snarls} -k {input.pack} -Aaz $RPATHS -s {wildcards.sample} {input.gbz} | gzip > {output}
        """

rule dv_make_examples:
    input:
        ref="results/pg/{graph}.ref.fa",
        ref_idx="results/pg/{graph}.ref.fa.fai",
        bam="results/bam/{sample}.{graph}.bam",
        bai="results/bam/{sample}.{graph}.bam.bai"
    output: "results/dv/make_examples.{sample}.{graph}.tfrecord.tar.gz",
    params:
        label="{sample}.{graph}",
        dvdir="dv_{sample}_{graph}"
    threads: 8
    container: "docker://google/deepvariant:1.5.0"
    benchmark: 'benchmark/results/dv/dv_make_examples.{sample}.{graph}.benchmark.tsv'
    shell:
        """
        mkdir -p {params.dvdir}

        seq 0 $(({threads}-1)) | \
        parallel -q --halt 2 --line-buffer /opt/deepvariant/bin/make_examples \
        --mode calling \
        --ref {input.ref} \
        --reads {input.bam} \
        --examples ./{params.dvdir}/make_examples.{params.label}.tfrecord@{threads}.gz \
        --sample_name {wildcards.sample} \
        --gvcf ./{params.dvdir}/gvcf.{params.label}.tfrecord@{threads}.gz \
        --channels insert_size \
        --min_mapping_quality 1 \
        --keep_legacy_allele_counter_behavior \
        --task {{}}
        
        tar -czvf '{output}' {params.dvdir}
        rm -rf {params.dvdir}
        """

rule dv_call_variants:
    input:
        ref="results/pg/{graph}.ref.fa",
        ref_idx="results/pg/{graph}.ref.fa.fai",
        ex="results/dv/make_examples.{sample}.{graph}.tfrecord.tar.gz"
    output:
        vcf="results/vcf/{sample}.{graph}.snv_indels.vcf.gz",
        gvcf="results/vcf/{sample}.{graph}.snv_indels.g.vcf.gz"
    params:
        call_tf="temp.call_variants_output.{sample}.tfrecord.gz",
        label="{sample}.{graph}",
        dvdir="dv_{sample}_{graph}"
    container: "docker://google/deepvariant:1.5.0"
    threads: 8
    benchmark: 'benchmark/results/vcf/{sample}.{graph}.dv_call_variants.benchmark.tsv'
    shell:
        """
        tar -xzf {input.ex}

        /opt/deepvariant/bin/call_variants \
        --outfile {params.call_tf} \
        --examples "{params.dvdir}/make_examples.{params.label}.tfrecord@{threads}.gz" \
        --checkpoint /opt/models/wgs/model.ckpt

        /opt/deepvariant/bin/postprocess_variants \
        --ref {input.ref} \
        --infile {params.call_tf} \
        --nonvariant_site_tfrecord_path "{params.dvdir}/gvcf.{params.label}.tfrecord@{threads}.gz" \
        --outfile "{output.vcf}" \
        --gvcf_outfile "{output.gvcf}"

        rm -fr {params.call_tf} {params.dvdir}
        """
