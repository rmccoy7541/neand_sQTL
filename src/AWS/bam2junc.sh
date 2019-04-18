#!/bin/bash

leafCutterDir='' ## use only if you don't have scripts folder in your path
bamfile=$1
bedfile=`basename $1`.bed
juncfile=$2
reffile=$3

if [ ! -z "$leafCutterDir" ]
then
    samtools view -T $reffile $bamfile | python $leafCutterDir/scripts/filter_cs.py | $leafCutterDir/scripts/sam2bed.pl --use-RNA-strand - $bedfile
    $leafCutterDir/scripts/bed2junc.pl $bedfile $juncfile
else
    if ! which filter_cs.py>/dev/null
    then 
        echo "ERROR:"
        echo "Add 'scripts' forlder to your path or set leafCutterDir variable in $0"
        exit 1
    fi
    samtools view -T $reffile $bamfile | filter_cs.py | sam2bed.pl --use-RNA-strand - $bedfile
    bed2junc.pl $bedfile $juncfile
fi
rm $bedfile