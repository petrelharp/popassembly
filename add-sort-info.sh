#!/bin/bash

IBDS=$(find . -regex "./chrom[0-9_]*/chrom_[0-9_a-z]*.[0-9]*.beagle.ibd")

for IBD in $IBDS
do
    # sort them and include the sort keys
    cat $IBD | awk ' BEGIN{OFS="\t"} { x=$1<$3?$1:$3;y=$1<$3?$3:$1; $9=sprintf("%06d%08d",x,y); print} ' | sort -n -k 9 > ${IBD}.sorted
done

