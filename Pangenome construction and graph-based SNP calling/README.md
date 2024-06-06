## Pipeline for constrction pangenome graph using minigraph-cactu and variant calling based on pangenome graph



1. Pear pangenome graph  construction 

```
cactus-pangenome ./js ./Pyrus.pg.txt --mgCores 100 --mapCores 100 --consCores 100 --gbz clip full filter  --gfa clip full filter --vcf clip full filter --odgi clip full filter --chrom-vg clip full filter --chrom-og clip full filter --outDir ./ --outName Pyrus.pg.output --reference YLXB
```



2. Using 'filter graph' as index to perform mapping step

```bash
cd index
ls 
Pyrus-PG.d2.dist  Pyrus-PG.d2.gbz  Pyrus-PG.d2.min Pyrus-PG.snarls Pyrus-PG.xg
cd ..
vg giraffe -Z index/Pyrus-PG.d2.gbz -m index/Pyrus-PG.d2.min -d index/Pyrus-PG.d2.dist --progress --fastq-in data/pear1_1.QC.fq.gz --fastq-in data/pear1_2.QC.fq.gz -N pear1 -t 100 -R pear1 > pear1.gam  2> pear1.err
vg surject -x index/Pyrus-PG.d2.xg -m pear1.gam --threads 20 --sam-output > pear1.sam
#fix bam header
samtools view -h pear1.bam |sed -e 's/YLXB#0#//g;'|samtools sort --threads 10 -m 2g -O BAM > pear1.sorted.bam
```



3.Small varient calling

!Note: Please use deepvariant:1.5.0 if you don't have gpu in your server

```bash
singularity run --nv -B ./:/data docker://google/deepvariant:1.5.0-gpu /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=/data/ref.fa --reads=/data/pear1.sorted.bam --output_vcf=/data/pear1.vcf.gz --output_gvcf=/data/pear1.g.vcf.gz --intermediate_results_dir=/data/tmp --logging_dir=/data/log --make_examples_extra_args="min_mapping_quality=1,keep_legacy_allele_counter_behavior=true,normalize_reads=true" --num_shards=10
```



4.SV calling

```bash
vg pack -t 10 -x index/Pyrus-PG.d2.gbz -g pear1.gam -o pear1.pack
vg call -t 10 index/Pyrus-PG.d2.gbz -r index/Pyrus-PG.snarls -k pear1.pack -s pear1 -a > pear1.vcf
```

