#! /bin/bash 

# load inStrain environment
export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate inStrain
PYTHONPATH=''

# arguments
mapping=$1
genome=$2

genomeName=$(basename $genome .fa)
mappingFile=$mapping-historicalVirusRefs.sorted.bam
mappingName=$mapping-inStrain

fasta=/home/GLBRCORG/emcdaniel/EBPR/Flanking/ref_genomes/all_genomes/$genomeName.fa
genes=/home/GLBRCORG/emcdaniel/EBPR/Flanking/ref_genomes/annotations/$genomeName.genes.fna

# cd to mapping folder
cd /home/GLBRCORG/emcdaniel/EBPR/Flanking/ref_genomes/relative_abundance/mappingResults/historical_virus_genomes

# profile command
inStrain profile $mappingFile $fasta -o ../../../inStrain/$mappingName/$genomeName-$mappingName.IS -p 8 -g $genes -s ../../../inStrain/all-flanking-historical-genomes.stb