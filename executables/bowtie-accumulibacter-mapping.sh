#! /bin/bash

# arguments
sample=$1
r1=$2
r2=$3
samplename=$sample-UW1-UW3-mapping

cd /home/GLBRCORG/emcdaniel/EBPR/Flanking/accumulibacter/mappingResults

# mapping command
/opt/bifxapps/bowtie2-2.3.5.1/bowtie2 --threads 4 -x /home/GLBRCORG/emcdaniel/EBPR/Flanking/accumulibacter/refGenomes/UW1_UW3_bt2/UW1_UW3.fasta -1 $r1 -2 $r2 > $samplename.sam


# BAM, sort, index
/opt/bifxapps/samtools-1.9/bin/samtools view -S -b $samplename.sam >  $samplename-spRep.bam
/opt/bifxapps/samtools-1.9/bin/samtools sort $samplename.bam -o  $samplename.sorted.bam
/opt/bifxapps/samtools-1.9/bin/samtools index $samplename-spRep.sorted.bam $samplename.sorted.bam.bai
