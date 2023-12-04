This repository contains a [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) workflow for the [vg toolkit](https://github.com/vgteam/vg).
With it, you can index pangenome, map short reads to the graph and call small variants.

## Prepare pangenome indexes

Starting from a GFA, the `index_pangenome` rule creates the GBZ, distance and minimizer indexes.

It starts from an *unzipped* GFA file, specified as the *gfa* field of the config.
A graph name (*graph*) is also provided to label the output files.

For example, to run the test data:

```
gunzip -c testdata/mhc.gfa.gz > testdata/mhc.gfa
snakemake --config gfa=testdata/mhc.gfa graph=mhc-test -p index_pangenome --cores 2
```

## Analyzing samples with short-read sequencing

Input FASTQ files are taken from columns *fq1* and *fq2* of the TSV specified by *sample_tsv* in the config.
If *sample* if specified in the config, only this sample, or these samples if it's an array or white-space-separated list of samples, will be run.
Otherwise, all samples in the *sample_tsv* TSV file will be analyzed.

See [config/sample_info.tsv](config/sample_info.tsv) for an example.

### Map short reads to a pangenome

```
snakemake --config gbz=testdata/mhc.gbz hapl=testdata/mhc.hapl ref_fa=testdata/mhc.ref.fa graph=mhc-test sample_tsv=config/sample_info.tsv sample=samp1 -p map_short_reads --cores 2 -n
```

If needed, the pangenome will be indexed with the rule described above.

### Call small variants with DeepVariant

Here, the reads will be projected to a linear reference genome represented by paths in the pangenome.
The user must specify a text file listing those reference paths (one path per line). 
This paths list file is specified in the config with *ref_paths_list*.

```
snakemake --config gfa=testdata/mhc.gfa graph=mhc-test sample_tsv=config/sample_info.tsv sample=samp1 ref_paths_list=testdata/mhc.paths_list.txt -p call_small_variants_from_short_reads --cores 2 --use-singularity -n

snakemake --config gbz=testdata/mhc.gbz hapl=testdata/mhc.hapl ref_fa=testdata/mhc.ref.fa ref_fa_idx=testdata/mhc.ref.fa.fai graph=mhc-test sample_tsv=config/sample_info.tsv sample=samp1 -p call_small_variants_from_short_reads --cores 2 --use-singularity -n
```

If needed, the pangenome will be indexed and the reads mapped with the rules described above.

### Genotype variants in the pangenome from short reads mapping

Same as above, the genotyped variants will ultimately be represented relative to a linear reference genome.
The user can specify the list of reference paths in the config with *ref_paths_list*, a text file listing the names of those paths.

```
snakemake --config gfa=testdata/mhc.gfa graph=mhc-test sample_tsv=config/sample_info.tsv sample=samp1 ref_paths_list=testdata/mhc.paths_list.txt -p genotype_variants_from_short_reads --cores 2 -n
```

If needed, the pangenome will be indexed and the reads mapped with the rules described above.

## Using existing pangenome indexes

Make sure the files are dated consistently with their dependencies. 
For example, if the GBZ file is more recent than the other indexes, the workflow will re-create them.
The order, for the oldest to most recent file is: `.gbz`, `.dist`, `.hapl`, `.ref.fa`, `.ref.fa.fai`.

### Human Pangenome Reference Consortium pangenomes

Download the Minigraph-Cactus from the [human-pangenomics/hpp_pangenome_resources repo](https://github.com/human-pangenomics/hpp_pangenome_resources#minigraph-cactus):

```sh
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.gbz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.hapl
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.snarls
```

We'll project the reads/variants on GRCh38 so we prepare a list of paths corresponding to the main chromosomes:

```sh
vg paths -RL -x hprc-v1.1-mc-grch38.gbz | grep GRCh38 | grep -v "_" | grep -v EBV > hprc-v1.1-mc-grch38.paths_list.txt
```

Run the workflow:

```sh
snakemake --configfile config/config.hprc.yaml -p genotype_variants_from_short_reads --cores 2 -n --slurm --use-singularity --profile profile/default
```
