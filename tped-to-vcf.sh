#!/bin/bash

# from ../../../thedata/POPRES_Genotypes_QC2_v2_TXT.tfam.gz
TFAM="../everyone.tfam"

if ! [ $# -eq 1 ]
then
    echo "Usage: $0 (.tped file)"
    exit 1
fi

TPED_FILE="$1"

if ! [ -f $TPED_FILE -a -f $TFAM ]
then
    echo "Can't find" $TPED_FILE
    exit 1
fi

TPED="${TPED_FILE%.tped.gz}"
# this makes a bunch of things we don't need; put them in here
TEMPDIR="${TPED}_dir"
mkdir -p $TEMPDIR

gzip -d $TPED_FILE
ln -s $TFAM ${TPED}.tfam
plink --tfile $TPED --noweb --make-bed --out $TEMPDIR/$TPED
gzip ${TPED_FILE%.gz}

cd $TEMPDIR || ( echo "can't change directories to $TEMPDIR"; exit 1 )
pseq $TPED new-project
pseq ${TPED}.pseq load-plink --file $TPED --id $TPED
pseq ${TPED}.pseq write-vcf | awk '$5 == "0" { $5="A" }; {print;}' | gzip -c >../${TPED}.vcf.gz
cd ..
rm -rf $TEMPDIR

echo "All done with $TPED"
exit 0
