#! /bin/bash

# Map metagenomic reads to refs/assemblies

export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate coverM
PYTHONPATH=''

BBPATH=/opt/bifxapps/bbmap-38.32
SAMTOOLS=/opt/bifxapps/samtools-1.9/bin

# cd to directory of metagenomes
cd /home/GLBRCORG/emcdaniel/EBPR/Flanking/metagenomes

ref=$1
meta=$2
out=$3

$BBPATH/bbmap.sh ref=$ref in=$meta outm=$out idtag minid=0.95 nodisk -Xmx50g
$SAMTOOLS/samtools sort $out -o ${out%.bam}.sorted.bam
