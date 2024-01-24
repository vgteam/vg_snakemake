rule pack:
    input:
        gbz="results/{sample}/{graph}.sample_pg.{sample}.gbz",
        gaf="results/{sample}/{sample}.{graph}.gaf.gz"
    output: temp('results/{sample}/{sample}.{graph}.pack')
    threads: 8
    priority: 3
    benchmark: 'benchmark/{sample}.{graph}.pack.benchmark.tsv'
    log: 'logs/pack.{graph}.{sample}.log'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    shell: "vg pack -x {input.gbz} -a {input.gaf} -Q 5 -t {threads} -o {output} 2> {log}"

rule vgcall:
    input: 
        gbz="results/{sample}/{graph}.sample_pg.{sample}.gbz",
        paths_list=config['ref_paths_list'],
        pack='results/{sample}/{sample}.{graph}.pack'
    output: 'results/{sample}/{sample}.{graph}.gt.minlen{minlen}.vcf.gz'
    threads: 8
    priority: 4
    benchmark: 'benchmark/{sample}.{graph}.vgcall.minlen{minlen}.benchmark.tsv'
    container: "docker://quay.io/vgteam/vg:v1.52.0"
    shell:
        """
        RPATHS=""
        for RP in `cat {input.paths_list}`
        do
        RPATHS=`echo "-p $RP $RPATHS"`
        done
        
        vg call -t {threads} -k {input.pack} -Az $RPATHS -s {wildcards.sample} -c {wildcards.minlen} {input.gbz} | gzip > {output}
        """

rule dv_make_examples:
    input:
        ref=getref(),
        ref_idx=getrefidx(),
        bam="results/{sample}/{sample}.{graph}.{surj}.bam",
        bai="results/{sample}/{sample}.{graph}.{surj}.bam.bai"
    output: temp("results/{sample}/{sample}.{graph}.{surj}.make_examples.tfrecord.tar.gz")
    params:
        label="{sample}.{graph}.{surj}",
        dvdir="dv_{sample}_{graph}_{surj}"
    threads: 8
    priority: 6
    container: "docker://google/deepvariant:1.5.0"
    benchmark: 'benchmark/{sample}.{graph}.{surj}.dv_make_examples.benchmark.tsv'
    log: 'logs/{sample}.{graph}.{surj}.dv_make_examples.log'
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
        --task {{}} 2> {log}
        
        tar -czf '{output}' {params.dvdir}
        rm -rf {params.dvdir}
        """

if config['use_gpu']:
    rule dv_call_variants_gpu:
        input:
            ref=getref(),
            ref_idx=getrefidx(),
            ex="results/{sample}/{sample}.{graph}.{surj}.make_examples.tfrecord.tar.gz"
        output:
            vcf="results/{sample}/{sample}.{graph}.{surj}.snv_indels.vcf.gz",
            gvcf="results/{sample}/{sample}.{graph}.{surj}.snv_indels.g.vcf.gz"
        params:
            call_tf="temp.call_variants_output.{sample}.{surj}.tfrecord.gz",
            label="{sample}.{graph}.{surj}",
            dvdir="dv_{sample}_{graph}_{surj}"
        container: "docker://google/deepvariant:1.5.0-gpu"
        threads: 8
        priority: 7
        benchmark: 'benchmark/{sample}.{graph}.{surj}.dv_call_variants_gpu.benchmark.tsv'
        log: 'logs/{sample}.{graph}.{surj}.dv_call_variants_gpu.log'
        shell:
            """
            tar -xzf {input.ex}

            /opt/deepvariant/bin/call_variants \
            --outfile {params.call_tf} \
            --examples "{params.dvdir}/make_examples.{params.label}.tfrecord@{threads}.gz" \
            --checkpoint /opt/models/wgs/model.ckpt 2> {log}

            /opt/deepvariant/bin/postprocess_variants \
            --ref {input.ref} \
            --infile {params.call_tf} \
            --nonvariant_site_tfrecord_path "{params.dvdir}/gvcf.{params.label}.tfrecord@{threads}.gz" \
            --outfile "{output.vcf}" \
            --gvcf_outfile "{output.gvcf}" 2>> {log}

            rm -fr {params.call_tf} {params.dvdir}
            """
else:
    rule dv_call_variants:
        input:
            ref=getref(),
            ref_idx=getrefidx(),
            ex="results/{sample}/{sample}.{graph}.{surj}.make_examples.tfrecord.tar.gz"
        output:
            vcf="results/{sample}/{sample}.{graph}.{surj}.snv_indels.vcf.gz",
            gvcf="results/{sample}/{sample}.{graph}.{surj}.snv_indels.g.vcf.gz"
        params:
            call_tf="temp.call_variants_output.{sample}.{surj}.tfrecord.gz",
            label="{sample}.{graph}.{surj}",
            dvdir="dv_{sample}_{graph}_{surj}"
        container: "docker://google/deepvariant:1.5.0"
        threads: 8
        priority: 7
        benchmark: 'benchmark/{sample}.{graph}.{surj}.dv_call_variants.benchmark.tsv'
        log: 'logs/{sample}.{graph}.{surj}.dv_call_variants.log'
        shell:
            """
            tar -xzf {input.ex}

            /opt/deepvariant/bin/call_variants \
            --outfile {params.call_tf} \
            --examples "{params.dvdir}/make_examples.{params.label}.tfrecord@{threads}.gz" \
            --checkpoint /opt/models/wgs/model.ckpt 2> {log}

            /opt/deepvariant/bin/postprocess_variants \
            --ref {input.ref} \
            --infile {params.call_tf} \
            --nonvariant_site_tfrecord_path "{params.dvdir}/gvcf.{params.label}.tfrecord@{threads}.gz" \
            --outfile "{output.vcf}" \
            --gvcf_outfile "{output.gvcf}" 2>> {log}

            rm -fr {params.call_tf} {params.dvdir}
            """

