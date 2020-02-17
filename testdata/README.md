## Test on a small dataset

We use the rules from the file in the repo's root:

```
ln -s ./Snakefile .
```

### Splitting reads in chunks

``` 
snakemake --configfile config.testdata.yaml -p genotype

## Diagrams of the jobs/rules/files
snakemake --configfile config.testdata.yaml --dag genotype | dot -Tsvg > ../imgs/construct-map-geno-dag.svg
snakemake --configfile config.testdata.yaml --filegraph genotype | dot -Tsvg > ../imgs/construct-map-geno-filegraph.svg
snakemake --configfile config.testdata.yaml --rulegraph genotype | dot -Tsvg > ../imgs/construct-map-geno-rulegraph.svg

## Clean up
snakemake --configfile config.testdata.yaml -p cleanall
```

### Mapping all reads in one job

Some times we might want to use only one job to map all the reads instead of splitting them in chunks.
To do that, set *nb_split_reads=0*.

``` 
snakemake --configfile config.testdata.yaml -p genotype --config nb_split_reads=0

## Clean up
snakemake --configfile config.testdata.yaml -p cleanall
```

### Giraffe mapper

The new and faster mapper uses different indexes.

``` 
snakemake --configfile config.testdata.yaml -p genotype

## Diagrams of the jobs/rules/files
snakemake --configfile config.testdata.yaml --dag genotype --config mapper="gaffe" | dot -Tsvg > ../imgs/construct-map-geno-dag.svg
snakemake --configfile config.testdata.yaml --filegraph genotype --config mapper="gaffe" | dot -Tsvg > ../imgs/construct-map-geno-filegraph.svg
snakemake --configfile config.testdata.yaml --rulegraph genotype --config mapper="gaffe" | dot -Tsvg > ../imgs/construct-map-geno-rulegraph.svg

## Clean up
snakemake --configfile config.testdata.yaml -p cleanall --config mapper="gaffe" 
```
