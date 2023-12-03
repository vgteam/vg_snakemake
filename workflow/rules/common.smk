import pandas as pd

if 'sample' not in config:
    config['sample'] = []

info = {}

if 'sample_tsv' in config:
    info = pd.read_csv(config["sample_tsv"], sep="\t", dtype={"sample": str}).set_index("sample", drop=False).sort_index()
    
    if len(config['sample']) == 0:
        config['sample'] = list(pd['sample'])

def getfq1(wildcards):
    # complain if sample not in the file
    if wildcards.sample not in info.fq1:
        print("Error: " + wildcards.sample + ' not in ' + config['sample_tsv'])
    # return fastq path for sample
    return info.fq1[wildcards.sample]

def getfq2(wildcards):
    # complain if sample not in the file
    if wildcards.sample not in info.fq2:
        print("Error: " + wildcards.sample + ' not in ' + config['sample_tsv'])
    # return fastq path for sample
    return info.fq2[wildcards.sample]

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

def getsnarls():
    if 'snarls' in config:
       return config['snarls']
    else:
        return "results/pg/{graph}.snarls"
