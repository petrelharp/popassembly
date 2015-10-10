#!/bin/bash
#PBS -S /bin/bash
#PBS -q cmb
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -l pmem=4gb
#PBS -l mem=4gb
#PBS -l vmem=4gb


if [ -e /usr/usc/java/1.8.0_45/setup.sh ]
then
    # on the cluster
    source /usr/usc/java/1.8.0_45/setup.sh 
    JAVA="java -Xmx3500m"
    BEAGLE="/home/rcf-40/pralph/cmb/software/beagle/beagle.09Oct15.56b.jar"
    export _JAVA_OPTIONS="-Xmx3500m -XX:CompressedClassSpaceSize=3500m"
else
    # at home
    JAVA="/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java"
    BEAGLE="/home/peter/software/beagle/beagle.09Oct15.56b.jar"

    if ! [ $# -eq 1 ]
    then
        echo "Usage: beagle.sh (.vcf file)"
        echo "or: qsub -vVCF_FILE='chrom_17_204992/chrom_204992_contig_15.vcf.gz' beagle.sh"
        exit 1
    fi

    VCF_FILE="$1"
fi

echo "VCF  file: $VCF_FILE"
echo "beagle: $BEAGLE"

if [[ -z ${VCF_FILE:-} ]]
then
    echo "Can't find $VCF_FILE or $BEAGLE"
    exit 1
fi

$JAVA -jar $BEAGLE gt=$VCF_FILE ibd=true out=${VCF_FILE%%.vcf.gz}.beagle
