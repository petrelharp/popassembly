#!/bin/bash

if [ -e /usr/usc/java/1.8.0_45/setup.sh ]
then
    # on the cluster
    source /usr/usc/java/1.8.0_45/setup.sh 
    JAVA=java
    BEAGLE="/home/rcf-40/pralph/cmb/software/beagle/beagle.09Oct15.56b.jar"
else
    # at home
    JAVA="/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java"
    BEAGLE="/home/peter/software/beagle/beagle.09Oct15.56b.jar"

    if ! [ $# -eq 1 ]
    then
        echo "Usage: beagle.sh (.vcf file)"
        exit 1
    fi

    VCF_FILE="$1"
fi

echo "VCF  file: $VCF_FILE"
echo "beagle: $BEAGLE"
$JAVA -jar $BEAGLE gt=$VCF_FILE ibd=true out=${VCF_FILE%%.vcf.gz}.beagle
