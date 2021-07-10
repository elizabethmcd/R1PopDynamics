#! /bin/bash 

# load inStrain environment
export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate inStrain
PYTHONPATH=''

# arguments
bam=$1
outfolder=$(basename $bam .sorted.bam)-quick-profile

# cd to mapping results folder
cd /home/GLBRCORG/emcdaniel/EBPR/Flanking/ref_genomes/relative_abundance/mappingResults/historical_virus_genomes


# inStrain quick profile command

inStrain quick_profile -p 2 -s ../../../inStrain/all-flanking-historical-genomes.stb -o ../../../profiles/$outfolder $bam ../../all-flanking-virus-historical-genomes.fasta