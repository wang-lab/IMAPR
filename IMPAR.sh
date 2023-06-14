#!/bin/bash

input_file="$1"

if [ -z "$input_file" ]; then
  echo "Error: No input file provided."
  exit 1
fi

sample_name=$(grep "^sample_name" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
input_format=$(grep "^input_format" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
tumor_input=$(grep "^tumor_input" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
normal_input=$(grep "^normal_input" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
out_prefix=$(grep "^out_prefix" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')

thread=$(grep "^thread" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
ram=$(grep "^ram" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')

gatk=$(grep "^gatk" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
picard=$(grep "^picard" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
samtools=$(grep "^samtools" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
bcftools=$(grep "^bcftools" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
hisat2=$(grep "^hisat2" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')

fasta_ref=$(grep "^fasta_ref" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
gtf_ref=$(grep "^gtf_ref" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
genelist_ref=$(grep "^genelist_ref" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
dbsnp_ref=$(grep "^dbsnp_ref" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
germline_ref=$(grep "^germline_ref" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
PON_ref=$(grep "^PON_ref" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
hisat_ref=$(grep "^hisat_ref" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
igg_ref=$(grep "^igg_ref" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
hla_ref=$(grep "^hla_ref" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
pseudo_ref=$(grep "^pseudo_ref" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
tcga_PON_ref=$(grep "^tcga_PON_ref" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
radar_ref=$(grep "^radar_ref" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
darned_ref=$(grep "^darned_ref" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
REDI_ref=$(grep "^REDI_ref" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')

if [ -z "$sample_name" ] || [ -z "$input_format" ] || [ -z "$tumor_input" ] || [ -z "$normal_input" ] || [ -z "$out_prefix" ] || [ -z "$thread" ] || [ -z "$ram" ] || [ -z "$gatk" ] || [ -z "$picard" ] || [ -z "$samtools" ] || [ -z "$bcftools" ] || [ -z "$hisat2" ] || [ -z "$fasta_ref" ] || [ -z "$gtf_ref" ] || [ -z "$genelist_ref" ] || [ -z "$dbsnp_ref" ] || [ -z "$germline_ref" ] || [ -z "$PON_ref" ] || [ -z "$hisat_ref" ] || [ -z "$igg_ref" ] || [ -z "$hla_ref" ] || [ -z "$pseudo_ref" ] || [ -z "$tcga_PON_ref" ] || [ -z "$radar_ref" ] || [ -z "$darned_ref" ] || [ -z "$REDI_ref" ]; then
  echo "Error: Found missing parameters in $input_file file."
  exit 1
fi

if ! command -v $gatk > /dev/null; then
  echo "Error: GATK took kit is not installed or not in the PATH."
  exit 1
fi

if ! command -v java -jar $picard > /dev/null; then
  echo "Error: Picard tool is not installed or not in the PATH."
  exit 1
fi

if ! command -v $samtools > /dev/null; then
  echo "Error: samtools tool is not installed or not in the PATH."
  exit 1
fi

if ! command -v $bcftools > /dev/null; then
  echo "Error: bcftools tool is not installed or not in the PATH."
  exit 1
fi

if ! command -v $hisat2 > /dev/null; then
  echo "Error: hisat2 is not installed or not in the PATH."
  exit 1
fi

echo "perl detect_variants.pl -ID $sample_name -mode $input_format -T $tumor_input -N $normal_input -R $fasta_ref -O $out_prefix -thread $thread -ram $ram -gatk $gatk -picard $picard -hisat2 $hisat2 -samtools $samtools -gtf $gtf_ref -gene $genelist_ref -dbsnp $dbsnp_ref -hisat2_reference $hisat_ref -germline $germline_ref -pon $PON_ref"
perl detect_variants.pl -ID $sample_name -mode $input_format -T $tumor_input -N $normal_input -R $fasta_ref -O $out_prefix -thread $thread -ram $ram -gatk $gatk -picard $picard -hisat2 $hisat2 -samtools $samtools -gtf $gtf_ref -gene $genelist_ref -dbsnp $dbsnp_ref -hisat2_reference $hisat_ref -germline $germline_ref -pon $PON_ref

echo "perl filter_variants.pl -ID $sample_name -O $out_prefix -R $fasta_ref -igg $igg_ref -hla $hla_ref -pseudo $pseudo_ref -tcga $tcga_PON_ref -radar $radar_ref -darned $darned_ref -redi $REDI_ref"
perl filter_variants.pl -ID $sample_name -O $out_prefix -R $fasta_ref -igg $igg_ref -hla $hla_ref -pseudo $pseudo_ref -tcga $tcga_PON_ref -radar $radar_ref -darned $darned_ref -redi $REDI_ref -samtools $samtools -bcftools $bcftools

echo "perl machine_learning.pl -ID $sample_name -O $out_prefix -gtf $gtf_ref"
perl machine_learning.pl -ID $sample_name -O $out_prefix -gtf $gtf_ref
