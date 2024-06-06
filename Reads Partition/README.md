## Pipeline to partition reads (Hi-C, RNA-seq, Genomics) 



This script was largely inspired by two nice pervious work for trio-bining assembly of two previous studies [yak and cattle genome](https://academic.oup.com/gigascience/article/9/4/giaa029/5815405), [cattle trio-bining genome](https://www.nature.com/articles/s41467-020-15848-y). In Cattle trio-bining genome, they used a _"A haplotype score for a pair was defined as the sum of the percent identity multiplied by match length for each read end (unmapped read ends were assigned a score of 0). "_  to partition reads. In yak and cattle genome, they used a similar way to partition reads, and the indicator to partition reads was "AS" tag. A read pairs first need to map to each haplotype and will get two alignment scores, if pairs with a higher score for one haplotype were considered breed specific and assigned to their respective haplotype. Pairs with a tied score were considered homozygous and assigned to both haplotypes for scaffolding.



### Main steps

1.align the reads to both haplotypes

2.run the script for both two bams generated from step1

3.Partitioned reads based on results of python script.



### Practice

1.mapping short reads to each maternal and paternal genomes/contigs 

```bash
bwa-mem2 index -p maternal maternal.genome.fa
bwa-mem2 index -p paternal paternal.genome.fa
bwa-mem2 mem maternal read_r1.fq.gz read_r2.fq.gz -t 10 | samtools view -bh - | samtools sort -n - > maternal.sorted.bam
bwa-mem2 mem paternal read_r1.fq.gz read_r2.fq.gz -t 10 | samtools view -bh - | samtools sort -n - > paternal.sorted.bam
```

2. Running the script for both two bams generated from step1

```bash
Classified_by_Alignment.py maternal.sorted.bam paternal.sorted.bam
```

3. Getting reads from raw fq file

```bash
seqtk subseq read_r1.fq.gz maternal.sorted.bam.id > maternal.read_r1.fq
seqtk subseq read_r2.fq.gz maternal.sorted.bam.id >maternal.read_r2.fq
#
seqtk subseq read_r1.fq.gz paternal.sorted.bam.id > paternal.read_r1.fq
seqtk subseq read_r2.fq.gz paternal.sorted.bam.id >paternal.read_r2.fq
```

