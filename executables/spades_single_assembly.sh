#! /bin/bash

######################
# Queue single metagenomic assemblies with SPAdes from list
######################

# Metagenomic reads path queued from provided file
# Make sure reads have been filtered previously before assembly - garbage in = garbage out
sample=$1
r1=$2
r2=$3
outdir=/home/GLBRCORG/emcdaniel/EBPR/Flanking/assemblies/new_assemblies

# Python path
#PYTHONPATH=/opt/bifxapps/bin/python3.4
    # use specific 3.4 python path for running SPADES
SPADESPATH=/opt/bifxapps/SPAdes-3.9.0-Linux/bin/

# Run spades
/opt/bifxapps/bin/python3.4 $SPADESPATH/spades.py -t 20 -m 800 -k 21,33,55,77,99,127 -1 $r1 -2 $r2 -o $outdir/$sample