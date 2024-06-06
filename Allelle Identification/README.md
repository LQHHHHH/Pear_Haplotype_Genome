## Pipeline to identify alleles from two haplotype genomes





####  1. Identify synteny genes from two genomes ---- gene level identification

```
#prepare Cds and bed
python -m jcvi.formats.gff bed --type=mRNA --key=ID hap1.gff -o hap1.bed
python -m jcvi.formats.gff bed --type=mRNA --key=ID hap1.gff -o hap2.bed
python -m jcvi.compara.catalog ortholog hap1 hap2 --no_strip_names --cscore=0.99 --full
```



#### 2. Idenity Syteny region via anchorwave. Running Anchorwave

```bash
ref_prefix=$1
qry_prefix=$2
ref=$3
qry=$4
refcds=$5
qrycds=$6
refgff=$7
self_sam="$ref_prefix".gmap.self.sam
compar_sam="$qry_prefix".vs."$ref_prefix".gmap.sam
gmap_build --dir=./ --genomedb=$ref_prefix $ref
gmap_build --dir=./ --genomedb=$qry_prefix $qry
echo "gmap -t 20 -A -f samse -d $ref_prefix -D $ref_prefix/ $refcds > $self_sam"
gmap -t 20 -A -f samse -d $ref_prefix -D $ref_prefix/ $refcds > $self_sam
echo "gmap -t 20 -A -f samse -d $qry_prefix -D $qry_prefix/ $refcds > $compar_sam"
gmap -t 20 -A -f samse -d $qry_prefix -D $qry_prefix/ $refcds > $compar_sam

echo "anchorwave proali -t 10 -i $refgff  -as $refcds -r $ref  -a $compar_sam -ar $self_sam -s $qry -R 2 -Q 2 -n "$ref_prefix"_vs_"$qry_prefix".gmap.anchor -o "$ref_prefix"_vs_"$qry_prefix".gmap.maf -f "$ref_prefix"_vs_"$qry_prefix".gmap.f.maf -m 0"
anchorwave proali -t 10 -i $refgff  -as $refcds -r $ref  -a $compar_sam -ar $self_sam -s $qry -R 2 -Q 2 -n "$ref_prefix"_vs_"$qry_prefix".gmap.anchor -o "$ref_prefix"_vs_"$qry_prefix".gmap.maf -f "$ref_prefix"_vs_"$qry_prefix".gmap.f.maf -m 0
```



#### 3.Run diamond for calculate sequence similarity

```bash
$diamond makedb --in hap1.pep -d hap1
$diamond blastp -d hap1 -q hap2.pep -o hap1_vs_hap2.ultrasens.blast -p 20 --evalue 1e-5
```



#### 4. Transform Anchorwave output to used format

```bash
$anchor.anchor2coord.py HXSA_vs_HXSB.gmap.anchor HXSA_vs_HXSB.gmap.anchor.coords
$sort -k10,10 -k11,11 -k1,1n -k2,2n -k3,3n -k4,4n HXSA_vs_HXSB.gmap.maf.coords > hap1_vs_hap2.gmap.maf.coords.sorted
```



#### 5.  Transferring Anchorwave coords file to sqlite3 database and blast file

```bash
$sqlite3 hap1_vs_hap2.maf.db
create table dat(qstart INTEGER, qend INTEGER, sstart INTEGER, send INTEGER, qlens INTEGER, slens INTEGER,identity INTEGER, qstand text , sstand text , qchr text, schr text);
.separator "\t"
.import hap1_vs_hap2.gmap.maf.coords.sorted
.quit
```



#### 6.Transferring blast file to sqlite3 database

```bash
$sqlite3 hap1_vs_hap2.ultrasens.blast.db
create table dat(g1 text, g2 text, simarity INTEGER, length INTEGER, mismatch INTEGER, gapopen INTEGER,qstart INTEGER, qend INTEGER , sstart INTEGER , send INTEGER, evalue INTEGER, bitscore INTEGER);
.separator "\t"
.import hap1_vs_hap2.ultrasens.blast dat
.quit
```



#### 7. Identification potential alleles

```bash
$python3 Identify_alleles.anchorwave.v2.3.py hap1.hap2.1x1.lifted.anchors coords.maf.db HXS-A_HXS-B.bed hap1_vs_hap2.ultrasens.blast.db hap1-hap2.alleles.Anchorwave.txt
```





Finally, You need to remove potential tandem duplicated genes and only keep to one-to-one
