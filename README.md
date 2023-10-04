# Revisiting genomes of non-model species with long reads yields new insights into their biology and evolution

## Nanopore reads

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
### Assemblies 

```sh
flye -o flye_out --nano-hq ont_reads.fastq.gz
flye -o flye_out --nano-corr ont_reads.fastq.gz
flye -o flye_out --keep-haplotypes --nano-hq ont_reads.fastq.gz
flye -o flye_out --keep-haplotypes --nano-corr ont_reads.fastq.gz
```
```sh

```
## PacBio HiFi reads

### *k*-mer analyses

```sh
mkdir tmp_smudge
ls hifi_reads.fastq.gz  > FILES
kmc -k27 -ci1 -cs10000 @FILES kmcdb tmp_smudge
kmc_tools transform kmcdb histogram kmcdb_k27.hist -cx10000

L=$(smudgeplot.py cutoff kmcdb_k27.hist L)
U=$(smudgeplot.py cutoff kmcdb_k27.hist U)
echo $L $U

kmc_tools transform kmcdb -ci"$L" -cx"$U" dump -s kmcdb_L"$L"_U"$U".dump
smudgeplot.py hetkmers -o kmcdb_L"$L"_U"$U" < kmcdb_L"$L"_U"$U".dump

smudgeplot.py plot kmcdb_L"$L"_U"$U"_coverages.tsv
```

### Assemblies

```sh
flye -o flye_out --pacbio-hifi hifi_reads.fastq.gz
flye -o flye_out --keep-haplotypes --pacbio-hifi hifi_reads.fastq.gz
```
``
