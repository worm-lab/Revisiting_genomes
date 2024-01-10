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

### Quality and length evaluation

```sh
NanoPlot --fastq ont_reads.fastq.gz -o nanoplot_ont --raw
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

### Quality and length evaluation

```sh
NanoPlot --fastq hifi_reads.fastq.gz -o nanoplot_hifi --raw
```

### Assemblies

```sh
flye -o flye_out --pacbio-hifi hifi_reads.fastq.gz
flye -o flye_out --keep-haplotypes --pacbio-hifi hifi_reads.fastq.gz
```

```sh
# for Romanomermis culicivorax
hifiasm -l 3 -o hifiasm_l3 hifi_reads.fastq.gz

# for Panagrolaimus sp. PS1159
hifiasm -l 3 --n-hap 3 -o hifiasm_l3_nhap3 hifi_reads.fastq.gz
hifiasm -l 0 --n-hap 3 -o hifiasm_l0_nhap3 hifi_reads.fastq.gz
```

## PacBio HiFi + Nanopore reads

```sh
# for Romanomermis culicivorax
hifiasm -l 3 --ul ont_reads.trimmed.min15kb.fastq.gz -o hifiasm_l3_ont15kb hifi_reads.fastq.gz

# for Panagrolaimus sp. PS1159
hifiasm -l 3 --n-hap 3 --ul ont_reads.trimmed.min30kb.fastq.gz -o hifiasm_l3_nhap3_ont30kb hifi_reads.fastq.gz
hifiasm -l 0 --n-hap 3 --ul ont_reads.trimmed.min30kb.fastq.gz -o hifiasm_l0_nhap3_ont30kb hifi_reads.fastq.gz
```

## Assembly post-processing

### Decontamination

```sh
blastn -query assembly.fasta -db nt -outfmt "6 qseqid staxids bitscore std sscinames scomnames" \
               -max_hsps 1 -evalue 1e-25 -out assembly.fasta.blast.out
```
```sh
# for Nanopore reads
minimap2 -ax map-ont assembly.fasta ont_reads.trimmed.fastq.gz | samtools sort > minimap2.assembly.bam

# for HiFi reads
minimap2 -ax map-hifi assembly.fasta hifi_reads.fastq.gz | samtools sort > minimap2.assembly.bam
```
```sh
docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v5.4.7_cv1 busco -i assembly.fasta -o busco_metazoa_odb10_assembly -m genome -l metazoa_odb10
```
```sh
blobtools add --fasta assembly.fasta --cov minimap2.assembly.bam --hits assembly.fasta.blast.out --busco busco_metazoa_odb10_assembly/run_metazoa_odb10/full_table.tsv --taxdump taxdump --create blobdir_out
blobtools view blobdir_out
```

### Haplotig purging

```sh
# for Nanopore reads
minimap2 -x map-ont assembly.fasta ont_reads.trimmed.fastq.gz | gzip -c - > minimap2.assembly.paf.gz

# for HiFi reads
minimap2 -x map-hifi assembly.fasta hifi_reads.fastq.gz | gzip -c - > minimap2.assembly.paf.gz
        
pbcstat minimap2.assembly.paf.gz 
calcuts PB.stat > cutoffs 2>calcults.log

hist_plot.py -c cutoffs PB.stat purge_dups.${pri_asm}.png

split_fa assembly.fasta > assembly.fasta.split
minimap2 -xasm5 -DP assembly.fasta.split assembly.fasta.split | gzip -c - > assembly.fasta.split.self.paf.gz

purge_dups -2 -T cutoffs -c PB.base.cov assembly.fasta.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs dups.bed assembly.fasta
```

## Assembly evaluation

```sh
docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v5.4.7_cv1 busco -i assembly.fasta -o busco_metazoa_odb10_assembly -m genome -l metazoa_odb10
docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v5.4.7_cv1 busco -i assembly.fasta -o busco_nematoda_odb10_assembly -m genome -l nematoda_odb10
```
```sh
meryl count k=27 hifi_reads.fastq.gz output hifi.meryl
merqury.sh hifi.meryl/ final_assembly.fasta merqury_out
```

## Repeat and gene prediction

### Repeat masking and transposable element annotation

```sh
EDTA.pl --genome final_assembly.fasta --sensitive 1 --anno 1 --overwrite 1 --force 1
```

### Gene prediction

```sh
braker.pl --species Species_name --genome final_assembly.fasta --gff3 --UTR off --bam hisat2.assembly.bam
```
