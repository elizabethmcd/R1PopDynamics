#! /bin/bash

###################
# Filtering metagenomic samples with fastp
# For use on WEI GLBRC servers running HT Condor
##################

# set path where fastp is installed in local home directory bin
FASTPATH=/home/GLBRCORG/emcdaniel/bin

# queueing r1 r2 metagenomic reads and output folder/file names
sample=$1
file1=$2
file2=$3
out1=$sample-R1.qced.fastq
out2=$sample-R2.qced.fastq

cd /home/GLBRCORG/emcdaniel/EBPR/Flanking/metagenomes/cleaned_fastqs

$FASTPATH/fastp --in1 $file1 --in2 $file2 --out1 $out1 --out2 $out2 --detect_adapter_for_pe  --cut_tail --cut_tail_window_size 15 --cut_tail_mean_quality 25 -h $sample.html