#!/bin/bash

TUMOR="/data3/gtang/git/data/tumor.bam"
NORMAL="/data3/gtang/git/data/normal.bam"
REF="/data1/gtang/genome/gdc/GRCh38.d1.vd1.fa"
OUTPRE="./out_put"

GATK="gatk"
PICARD="~/picard.jar"
HISAT2="hisat2"

GTF="/data1/gtang/genome/gdc/gencode.v36.annotation.gtf"
DBSNP="/data1/gtang/genome/gdc/common_all_20180418_update.vcf.gz"
HISAT2_REF="/data1/gtang/genome/gdc/hisat2_GRCh38"

if [ ! -e "$TUMOR" ]
then
    echo "File $TUMOR for -Tumor input does not exist. Please check your files."
    exit 1
fi

if [ ! -e "$NORMAL" ]
then
    echo "File $NORMAL for -Normal input does not exist. Please check your files."
    exit 1
fi

if [ ! -e "$REF" ]
then
    echo "File $REF for -R input does not exist. Please check your files."
    exit 1
fi

if [ ! -e "$GTF" ]
then
    echo "File $GTF for -gtf input does not exist. Please check your files."
    exit 1
fi

if [ ! -e "$DBSNP" ]
then
    echo "File $DBSNP for -dbsnp input does not exist. Please check your files."
    exit 1
fi

#Detect variants using RNA-seq data
perl p1.pl -Tumor $TUMOR -Normal $NORMAL -R $REF -outPrefix $OUTPRE -gatk $GATK -picard $PICARD -hisat2 $HISAT2 -gtf $GTF -dbsnp $DBSNP -hisat2_reference $HISAT2_REF

#Filter somatic mutations
perl p2.pl -inputDir $OUTPRE

#Machine learning improvement
perl p3.pl -inputDir $OUTPRE