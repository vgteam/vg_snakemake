import pandas as pd

if 'sample1' in config and config['sample'] is None:
    config['sample'] = []

# parse sample list: if specified and one item, try to split at white spaces
if 'sample' in config and not isinstance(config['sample'], list):
    config['sample'] = config['sample'].split(' ')

# init config to avoid errors
for cc in ['sample', 'refsynt_fa', 'adapters_fa', 'adapters_tsv',
           'gfa', 'seqn_prefix']:
    if cc not in config:
        config[cc] = []

# set which bam should be produced based on if indel realignment is on
# default is to use directly surjected reads
config['bam_mode'] = 'surj'
if 'indel_realign_reads' in config and config['indel_realign_reads']:
    # if indel realignment is on, use surjected reads that are then realigned
    config['bam_mode'] = 'surj_realn'

# should we try to use GPU? by default, no
if 'use_gpu' not in config:
    config['use_gpu'] = False

# should we try to use GPU? by default, no
if 'cram_ref' not in config:
    config['cram_ref'] = 'reference_for_CRAM.fasta'

if 'max_samps' in config:
    config['max_samps'] = int(config['max_samps'])

# load information about the samples (if a sample TSV is provided)
info = {}
if 'sample_tsv' in config:
    info = pd.read_csv(config["sample_tsv"], sep="\t", dtype={"sample": str}).set_index("sample", drop=False).sort_index()
    
    if len(config['sample']) == 0:
        config['sample'] = list(info['sample'])

# rules to either get the files specified by the user, or path that will trigger rules to make them
def getfq1(wildcards):
    # complain if sample not in the file
    if 'fq1' not in info or wildcards.sample not in info.fq1:
        # print("Error: " + wildcards.sample + ' not in ' + config['sample_tsv'])
        return 'results/{sample}/{sample}.1.fastq.gz'
    # return fastq path for sample
    return info.fq1[wildcards.sample]

def getfq2(wildcards):
    # complain if sample not in the file
    if 'fq2' not in info or wildcards.sample not in info.fq2:
        # print("Error: " + wildcards.sample + ' not in ' + config['sample_tsv'])
        return 'results/{sample}/{sample}.2.fastq.gz'
    # return fastq path for sample
    return info.fq2[wildcards.sample]

def getcram(wildcards):
    # complain if sample not in the file
    if 'cram' not in info or wildcards.sample not in info.cram:
        # print("Error: " + wildcards.sample + ' not in ' + config['sample_tsv'])
        return 'results/{sample}/{sample}.cram'
    # return fastq path for sample
    return info.cram[wildcards.sample]

def getgbz():
    if 'gbz' in config:
        return config['gbz']
    else:
        if 'gfa' not in config:
            print("Error: neither 'gbz' nor 'gfa' are specified.")
        else:
            return "results/pg/{graph}.gbz"

def gethapl():
    if 'hapl' in config:
        return config['hapl']
    else:
        return "results/pg/{graph}.hapl"

def getref():
    if 'ref_fa' in config:
       return config['ref_fa']
    else:
        return "results/pg/{graph}.ref.fa"

def getrefidx():
    if 'ref_fa' in config:
       return config['ref_fa'] + '.fai'
    else:
        return "results/pg/{graph}.ref.fa.fai"

def getrefdict():
    if 'ref_fa' in config:
        ref_fa = config['ref_fa']
        if ref_fa.endswith('.fa'):
            ref_dict = ref_fa[:-3] + '.dict'
        if ref_fa.endswith('.fasta'):
            ref_dict = ref_fa[:-6] + '.dict'
        return ref_dict
    else:
        return "results/pg/{graph}.ref.dict"

if 'rm_all_on_success' in config and config['rm_all_on_success']:
    def tempCond(filen):
        return (temp(filen))
else:
    def tempCond(filen):
        return (filen)
