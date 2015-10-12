#!/bin/bash

# Run this to find finished (stdout) or unfinished (stderr) runs.

VCFS=$(find . -regex "./chrom_[0-9_]*/chrom_[a-zA-Z0-9_]*.vcf.gz")

for VCF in $VCFS
do
    NAME=${VCF%.vcf.gz}
    FINISHED=$(grep -l "beagle.* finished" ${NAME}.*beagle.log)
    if [ -z "$FINISHED" ]
    then
        >&2 echo "$NAME :: no finished runs"
    else
        echo "$NAME :: "$FINISHED
    fi
done
