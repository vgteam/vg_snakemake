## Small pangenome for the MHC region

Built with Minigraph-Cactus following [those instructions](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md#mhc-graph).
We just decreased the number of haplotypes to 10 to speed up construction.

```sh
# download the sequences
wget -q https://zenodo.org/record/6617246/files/MHC-61.agc
# make the seqfile
mkdir -p mhc-fa ; rm -f mhc-seqfile.txt
for s in $(agc listset MHC-61.agc | head); do printf "${s}\tmhc-fa/${s}.fa\n" >> mhc-seqfile.txt; agc getset MHC-61.agc $s > mhc-fa/${s}.fa; done
# (Optional) clean up the sample names. if you don't do this, adjust --reference accordingly below
sed -i mhc-seqfile.txt -e 's/^MHC-00GRCh38/MHC-GRCh38/g' -e 's/^MHC-CHM13.0/MHC-CHM13/g'
# (Optional) clean up the contig names (replacing numbers with "MHC")
for f in mhc-fa/*.fa; do sed -i ${f} -e 's/#0/#MHC/g' -e 's/#1/#MHC/g' -e 's/#2/#MCH/g'; done
## run Minigraph-Cactus, eventually within docker container
## > docker run -it -v `pwd`:/app -w /app -u `id -u $USER` quay.io/comparative-genomics-toolkit/cactus:v2.6.13
cactus-pangenome ./js ./mhc-seqfile.txt --outDir mhc-pg --outName mhc --reference MHC-GRCh38 --mapCores 1
cp mhc-pg/mhc.gfa.gz mhc.gfa.gz
```

Prepare a file with the list of reference path names:

```sh
echo MHC-GRCh38 > mhc.paths_list.txt
```

Short reads were simulated from the pangenome with:

```sh
vg gbwt -o mhc.gbwt -g mhc.gg -Z mhc.gbz
vg sim -x mhc.gbz -l 150 -n 150000 -I -p 400 -v 20 -e .002 -i .00001 -m MHC-HG00438 -g mhc.gbwt -a | vg view -aX - > samp1.fastq
seqtk seq samp1.fastq -1 | awk '{gsub("_1$", ""); print $0}' | gzip > samp1_1.fastq.gz
seqtk seq samp1.fastq -2 | awk '{gsub("_2$", ""); print $0}' | gzip > samp1_2.fastq.gz
```

Test files: 

- `mhc.gfa.gz`: gzipped GFA file
- `mhc.gbz`: GBZ pangenome index 
- `mhc.min`: minimizer index for vg giraffe
- `mhc.dist`: distance index for vg giraffe
- `mhc.paths_list.txt`: text file listing the reference paths' names
- `samp1_1.fastq.gz`: first short read of the pair (simulated)
- `samp1_2.fastq.gz`: first short read of the pair (simulated)


#### CRAM file

Sometimes, we might want to use CRAM files as input.

```sh
samtools import -O CRAM -o samp1.cram samp1_1.fastq.gz samp1_2.fastq.gz
```
