# Revisiting genomes of non-model species with long reads yields new insights into their biology and evolution

## Nanopore assemblies

### Basecalling

```sh
pod5 convert fast5 fast5/*.fast5 --output pod5/
```
```sh
dorado duplex dna_r10.4.1_e8.2_400bps_sup\@v4.2.0 pod5/ > ont_reads.bam
```
```sh
samtools fastq ont_reads.bam > ont_reads.fastq.gz
```

### Trimming and read selection

```sh
# only trimming, no selection on quality
cat ont_reads.fastq | chopper > ont_reads.trimmed.q20.fastq

# select reads with a quality over 20
cat ont_reads.fastq | chopper -q 20 > ont_reads.trimmed.q20.fastq

gzip ont_reads*.fastq
```

