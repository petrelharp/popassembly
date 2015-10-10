#!/bin/bash

# from ../../../thedata/POPRES_Genotypes_QC2_v2_TXT.tfam.gz
TFAM="$(readlink -f $(dirname $0)/everyone.tfam)"

if ! [ $# -eq 1 ]
then
    echo "Usage: $0 (.tped file)"
    exit 1
fi

TPED_FILE="$1"

if ! [ -f $TPED_FILE -a -f ${TFAM:-} ]
then
    echo "Can't find $TPED_FILE or $TFAM"
    exit 1
fi

TPED="${TPED_FILE%.tped.gz}"
TPED_NAME=$(basename $TPED)
# this makes a bunch of things we don't need; put them in here
TEMPDIR="${TPED}_dir"
mkdir -p $TEMPDIR

gzip -d $TPED_FILE
ln -s $TFAM ${TPED}.tfam
plink --tfile $TPED --noweb --make-bed --out $TEMPDIR/$TPED_NAME
gzip ${TPED_FILE%.gz}

pushd $TEMPDIR || ( echo "can't change directories to $TEMPDIR"; exit 1 )
pseq $TPED_NAME new-project
pseq ${TPED_NAME}.pseq load-plink --file $TPED_NAME --id $TPED_NAME
pseq ${TPED_NAME}.pseq write-vcf | awk '$5 == "0" { $5="A" }; {print;}' | gzip -c >../${TPED_NAME}.vcf.gz
popd
rm -rf $TEMPDIR

echo "All done with $TPED"
exit 0
