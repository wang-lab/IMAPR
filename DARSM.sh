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

gatk=$(grep "^gatk" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
picard=$(grep "^picard" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
samtools=$(grep "^samtools" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
hisat2=$(grep "^hisat2" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')

fasta_ref=$(grep "^fasta_ref" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
gtf_ref=$(grep "^gtf_ref" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
genelist_ref=$(grep "^genelist_ref" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
RNA_edits_ref=$(grep "^RNA_edits_ref" $input_file | awk -F '\t' '{print $2}' | tr -d '\n')
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

if [ -z "$sample_name" ] || [ -z "$input_format" ] || [ -z "$tumor_input" ] || [ -z "$normal_input" ] || [ -z "$out_prefix" ] || [ -z "$gatk" ] || [ -z "$picard" ] || [ -z "$samtools" ] || [ -z "$hisat2" ] || [ -z "$fasta_ref" ] || [ -z "$gtf_ref" ] || [ -z "$genelist_ref" ] || [ -z "$RNA_edits_ref" ] || [ -z "$dbsnp_ref" ] || [ -z "$dbsnp_ref" ] || [ -z "$hisat_ref" ] || [ -z "$germline_ref" ] || [ -z "$PON_ref" ]; then
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

if ! command -v $hisat2 > /dev/null; then
  echo "Error: hisat2 is not installed or not in the PATH."
  exit 1
fi

echo "perl p1.pl -ID $sample_name -mode $input_format -T $tumor_input -N $normal_input -R $fasta_ref -O $out_prefix -gatk $gatk -picard $picard -hisat2 $hisat2 -gtf $gtf_ref -dbsnp $dbsnp_ref -hisat2_reference $hisat_ref"
#perl p1.pl -ID $sample_name -mode $input_format -T $tumor_input -N $normal_input -R $fasta_ref -O $out_prefix -gatk $gatk -picard $picard -hisat2 $hisat2 -gtf $gtf_ref -gene $genelist_ref -dbsnp $dbsnp_ref -hisat2_reference $hisat_ref -germline $germline_ref -pon $PON_ref

echo ""

echo "p2.pl -ID $sample_name -O $out_prefix -R $fasta_ref -igg $igg_ref -hla $hla_ref -pseudo $pseudo_ref -tcga $tcga_PON_ref -radar $radar_ref -darned $darned_ref -redi $REDI_ref"
perl p2.pl -ID $sample_name -O $out_prefix -R $fasta_ref -igg $igg_ref -hla $hla_ref -pseudo $pseudo_ref -tcga $tcga_PON_ref -radar $radar_ref -darned $darned_ref -redi $REDI_ref

