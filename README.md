This repository contains a [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) workflow for the [vg toolkit](https://github.com/vgteam/vg).
With it, you can: 

- index a pangenome with vg
- map short reads to the pangenome with vg Giraffe
- call small variants with [DeepVariant](https://github.com/google/deepvariant)
- genotype structural variants with vg call
- call structural variants with [manta](https://github.com/Illumina/manta)

If running the *all* rule, all the analysis listed above will be performed.
Specific analysis can also be performed using dedicated rules described below. 
Of note, the pipeline will automatically run the necessary steps for each analysis, which means, for example, that the read mapping will automatically be performed (if necessary) when using a variant calling rule.

Note: for the older version of the Snakemake pipeline see the [vg_genotyping_2020 branch](https://github.com/vgteam/vg_snakemake/tree/vg_genotyping_2020).

## Prepare pangenome indexes

This step is not necessary when using an existent pangenome (e.g. the HPRC pangenome, see below). 

Starting from a GFA, the `index_pangenome` rule creates the GBZ, distance and minimizer indexes.

It starts from an *unzipped* GFA file, specified as the *gfa* field of the config.
A graph name (*graph*) is also provided to label the output files.

For example, to run the test data:

```
gunzip -c testdata/mhc.gfa.gz > testdata/mhc.gfa
snakemake --config gfa=testdata/mhc.gfa graph=mhc-test -p index_pangenome --cores 2
```

The indexes will be placed in the `results/pg` directory.

## Analyzing samples with short-read sequencing

There are several ways to analyze short-read sequencing data.
The inputs might be in gzipped FASTQ files, but they might also be in CRAM files (e.g. to reanalyze data from existing projects).
Some steps are optional and can be enabled/disabled, for example the trimming of the raw reads, realignment of aligned reads, the use of GPU.

We'll start by describing how to specify the input files and the main default analysis.
At the end, we'll explain how to switch on/off specific options and run the workflow in HPCs.

### Input FASTQ/CRAM

The path to the input files should be compiled in a TSV file with a *sample* column, with the sample name, and other columns informing the path to the files for each sample.
This TSV file is specified by *sample_tsv* in the config.

Input FASTQ files are taken from columns *fq1* and *fq2* of the TSV specified by *sample_tsv* in the config.
See [config/sample_info.tsv](config/sample_info.tsv) and [config/config.cram.yaml](config/config.cram.yaml) for an example.

It is also possible to start from a CRAM file. 
In that case, a *cram* column should have the path to the CRAM file in the TSV file.
In addition, the reference FASTA used to make the CRAMs must be specified in the config with *cram_ref*.
See [config/sample_info.cram.tsv](config/sample_info.cram.tsv) and [config/config.cram.yaml](config/config.cram.yaml) for an example.
Note: *fq1*/*fq2* and *cram* columns can all be present in the TSV but, if both are available for a sample, the FASTQs will be used.

If *sample* if specified in the config, only this sample, or these samples (if it's an array or white-space-separated list of samples), will be run.
Otherwise, all samples in the *sample_tsv* TSV file will be analyzed.


### Map short reads to a pangenome

This will map short reads to the pangenome using vg Giraffe and produce a gzipped GAF file available at `results/{sample}.{graph}.gaf.gz`

```
snakemake --configfile config/config.yaml sample=samp1 -p map_short_reads --cores 2 -n
```

(remove `-n`, which means "dry-run", to actually run the workflow)

### Call small variants with DeepVariant

Here, the reads will be projected to a linear reference genome represented by paths in the pangenome.
The user can (should?) specify a text file listing those reference paths (one path per line). 
This paths list file is specified in the config with *ref_paths_list*.

```
snakemake --configfile config/config.yaml -p call_small_variants_from_short_reads --cores 2 --use-singularity -n
```

If needed, the reads will be mapped to the pangenome first.

### Genotype variants in the pangenome from short reads mapping

Same as above, the genotyped variants will ultimately be represented relative to a linear reference genome.
The user can specify the list of reference paths in the config with *ref_paths_list*, a text file listing the names of those paths.

```
snakemake --configfile config/config.yaml -p genotype_variants_from_short_reads --cores 2 -n
```

If needed, the pangenome will be indexed and the reads mapped with the rules described above.

Sometimes, we might want to genotype and report all variants sites, even if homozygous for the reference allele.
For example, to merge VCFs for multiple samples.
This mode can be enabled by setting `gt_ref` to `True` in the config (either in the `config.yaml` file, or `--config gt_ref=True` in the command line).
If `gt_ref` is *true*, `vg call` will use the `-a` parameter to force reporting on every snarl.

### Run all the analysis described above

```
snakemake --configfile config/config.yaml -p all --cores 2 -n
```

## Using existing pangenome indexes

In general, make sure the files are dated consistently with their dependencies. 
For example, if the GBZ file is more recent than the other indexes, the workflow will re-create them.
The order, from the oldest to most recent file, is: `.gbz`, `.dist`, `.hapl`, `.ref.fa`, `.ref.fa.fai`.

### Human Pangenome Reference Consortium pangenomes

Download the Minigraph-Cactus from the [human-pangenomics/hpp_pangenome_resources repo](https://github.com/human-pangenomics/hpp_pangenome_resources#minigraph-cactus):

```sh
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.gbz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.hapl
```

We'll project the reads/variants on GRCh38 so we've prepared a list of paths corresponding to the main chromosomes (see [config/hprc-v1.1-mc-grch38.paths_list.txt](config/hprc-v1.1-mc-grch38.paths_list.txt)).
This list was made using: 

```sh
vg paths -RL -x hprc-v1.1-mc-grch38.gbz | grep GRCh38 | grep -v "_" | grep -v EBV > config/hprc-v1.1-mc-grch38.paths_list.txt
```

The path names were also reordered (manually) because this order define the order in the BAMs.

By pointing to those file in the config (see [config/config.hprc.yaml](config/config.hprc.yaml)), the workflow can be run with:

```sh
snakemake --configfile config/config.hprc.yaml --cores 8 --use-singularity -p all -n
```

The `config/config.hprc.yaml` config file will of course need to be tweaked to specify your sample information (`sample_tsv` field).

## Understanding the config file

Main information to provide through the config (file or command line):

- `graph`: a name/label for the pangenome graph
- `gfa`: a pangenome file in GFA format. Must be unzipped.
- `gbz`: a pangenome file in GBZ format. At least `gfa` or `gbz` must be specified.
- `hapl`: a kmer index for each haplotype for the pangenome. Used by Giraffe to map read on a haplotype-sampled pangenome. Will be recomputed if missing.
- `ref_paths_list`: text file listing the paths in the pangenome to use as reference when projecting the reads or calling variants.
- `seqn_prefix`: the prefix to remove from these reference paths in the final BAM/VCF files. See *Chromosome re-naming* below.
- `sample_tsv`: the TSV file defining the sample and their input files (FASTQ, CRAM). See *Input FASTQ/CRAM* section above.
- `sample`: which samples to analyze? Either a sample name, a string with sample names separated by a white space, or an array of sample names. If missing or empty, all samples will be processed.

Other options:

- `ref_fa`: the linear reference genome sequence corresponding to the reference paths specified by`ref_paths_list`. Will be extracted from the pangenome if missing. See *Linear reference sequence and indexes* below.
- `gt_min_var_len`: the minimum length for variants to be genotyped with vg. Default is 30 (bp). Decrease, down to 0, to genotype small variants too.
- `indel_realign_reads`: should read be realigned to improve indel detection? See *Read re-alignment* section below.
- `cram_ref`: a reference for input CRAMs. See *Input FASTQ/CRAM* section above.
- `refsynt_fa`/`adapters_fa`/`adapters_tsv`: synthetic contaminants and adapter information to trim the raw FASTQ. See *Read trimming* section below.
- `use_gpu` should GPU be used for DeepVariant? See *Using GPUs* section below to set up GPUs.
- `rm_all_on_success`: should the results files be removed at the end? Useful to control disk space usage when paired with an archiving script. See *Minimizing local disk usage* section below.

### Linear reference sequence and indexes

The final steps of the pipeline, like variant calling, use a linear reference genome information. 
Specifically, they use a FASTA file, an FAI index for that file, and a sequence dictionnary.
If not provided by the user, the sequence will be extracted from the pangenome and indexed. 
Although this is fast, the user can still provide a reference FASTA file and indexes. 
The FASTA file is specified in `ref_fa` of the config. 

The index should be located in the same directory as the FASTA file and named appropriately, for example `ref.fa.fai` and `ref.dict` for a `ref.fa` FASTA file.

This is only recommended if those files are already available. 
If not, no need to prepare them, the pipeline will do it.

### Read trimming

If adapters and synthetic contaminant are provided, the raw reads in the FASTQs will first be trimmed.
QC of the FASTQ files before and after trimming is also performed.

The synthetic contaminant sequence(s) is provided as a FASTA file using *refsynt_fa* in the config.

The adapters sequence(s) is provided in both FASTA format and a TSV file with *adapters_fa* and *adapters_tsv*.

See examples of these files in [`config/config.trimming.yaml`](config/config.trimming.yaml), [`testdata/refsynt.fa`](testdata/refsynt.fa), [`testdata/adapters-test.fa`](testdata/adapters-test.fa), and [`testdata/adapters-test.tsv`](testdata/adapters-test.tsv).

### Read re-alignment

The reads can be realigned after projection to the linear reference genome. 
We use [ABRA2](https://github.com/mozack/abra2) for this step.
It improves the indel called by DeepVariant but increases the runtime significantly.
To switch the realignment off, set *indel_realign_reads* to *False* in the config.

It is only recommended if interested in the best performance for indel discovery.
If the main goal of the analysis is to look for SNVs, or genotype SVs, *indel_realign_reads* should be set to *False*.
Worst case, the BAMs could always be realigned later if really needed.

### Chromosome renaming

The *seqn_prefix* value, in the config, specifies the prefix that should be removed from the sequence names. 
For example, the HPRC pangenomes uses paths names that look like `GRCh38#0#chr1` to inform the genome of origin, haplotype, and contig name. 
Hence, in the config file (e.g. [`config/config.hprc.yaml`](config/config.hprc.yaml), we use *GRCh38#0#* as *seqn_prefix*. 

Caution: if provided as input files, the reference FASTA, index and dict, should already have this prefix stripped. 
The workflow won't strip this sequence prefix from user-provided files.

## Recommendations to adapt the workflow to an HPC

In practice, we can run one job (with very low resources) where we call the `snakemake` command and that will launch more jobs for each task.
The minimal dependencies would be Python (with snakemake and pandas installed) and Singularity.
Then, the `snakemake` command could look like:

```sh
snakemake --configfile my.config.yaml --slurm --profile profile/default --use-singularity -p all
```

In some cases, we might want to use tools already installed on the HPC, for example by loading "modules", instead of the singularity/docker containers.
This is not implemented yet, but could be added if there is interest.

## Using GPUs

If DeepVariant is run through containers with Singularity (default currently), a parameter needs to be added to enable GPU usage.
In addition to `--use-singularity`, we should add `--singularity-args '--nv -B .:/dum'` to the snakemake command.
This will make sure that the singularity instances are run with `--nv`.
(The `-B .:/dum` part is just a hack to avoid [a bug in snakemake](https://github.com/snakemake/snakemake/issues/1763).)

In addition, the value of `use_gpu` must be set to `True` in the config, either by editing the YAML config file passed with `--configfile`, or by adding `--config use_gpu=True` to the snakemake command.

In summary, the snakemake command, when using GPUs, looks like this:

```sh
snakemake --configfile my.config.yaml -p all --use-singularity --singularity-args '--nv -B .:/dum' --config use_gpu=True 
```

Using Slurm, it would look like this:

```sh
snakemake --configfile my.config.yaml -p all --slurm --profile profile/default --use-singularity --singularity-args '--nv -B .:/dum' --config use_gpu=True 
```

Activating GPUs on Slurm will depend on the slurm environment and vary from platform to platform.
To adapt to your Slurm environment, modify `dv_call_variants_gpu:slurm_extra` in the profile configuration file (passed with `--profile`, see [default config](profile/default/config.yaml)).
This parameter specifies which parameter should be passed to Slurm run a job with GPUs.
For example, it's currently `dv_call_variants_gpu:slurm_extra=--partition=gpuq --gres=gpu:A100_1g.10gb:1` because we need to include `--partition=gpuq --gres=gpu:A100_1g.10gb:1` to use special *gpuq* queue and specify GPUs in our Slurm environment ([Genotoul](https://bioinfo.genotoul.fr/)).

## Minimizing local disk usage

The first way to minimize disk usage is to make sure temporary files don't accumulate.
For this reason, temporary files have been marked so that they will be removed as soon as they are not needed. 

Secondly, if we have a large number of samples, we would prefer to run them in batches rather than starting all the samples at once.
Higher priority has been given to tasks later in the pipeline so jobs that continue/finish ongoing samples will be preferred over jobs starting new ones.
For the same reason , we can also force the pipeline to not analyze more than N samples at a time by creating fake dependencies. 
See an example in [`workflow/Snakefile_disk_opt_slurm`](workflow/Snakefile_disk_opt_slurm) and the *max_samps* config variable it uses.

Another approach might be to migrate/archive the results once they are ready (and delete the local copy from the server).
To do this, one could tweak the main `Snakefile`, adding a rule to perform this file transfer. *Example coming soon.*
For this mode, we can use `rm_all_on_success=True` in the Snakemake configuration to make sure even the "final" output files are removed as soon as they are not needed (in this scenario, once they've been migrated/archived).

## Test the workflow locally using the test data

```sh
snakemake --use-singularity --configfile config/config.yaml -p all --cores 2 -n
```
