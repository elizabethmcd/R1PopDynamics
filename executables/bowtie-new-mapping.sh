#! /bin/bash

# Map metagenomic reads to concatenated references with bowtie2 with newly assembled accumulibacter references

# Environment with samtools
export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate coverM
PYTHONPATH=''

# Arguments
sample=$1
samplename=$sample-newRefs

# cd to directory to output mapping results
cd /home/GLBRCORG/emcdaniel/EBPR/Flanking/ref_genomes/relative_abundance/mappingResults/new_genomes

# mapping command
/opt/bifxapps/bowtie2-2.3.5.1/bowtie2 --threads 6 -x /home/GLBRCORG/emcdaniel/EBPR/Flanking/ref_genomes/relative_abundance/new_bt2/new-bins.fasta --interleaved $sample > $samplename.sam


# BAM, sort, index
samtools view -S -b $samplename.sam >  $samplename.bam
samtools sort $samplename.bam -o  $samplename.sorted.bam
samtools index $samplename.sorted.bam $samplename.sorted.bam.bai