#!/bin/bash

input_folder="$1"
tools_folder="$2"

if [ -z "$input_folder" ]; then
  echo "Error: No output reference folder were provided, please specifiy where you want to store the reference file."
  exit 1
fi

if [ ! -d "$input_folder" ]; then
	mkdir $input_folder
fi

#fasta
echo "PREPARE fasta sequence file\n"
wget -O ./$input_folder/GRCh38.d1.vd1.fa.tar.gz https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834
tar -xzf ./$input_folder/GRCh38.d1.vd1.fa.tar.gz -C ./$input_folder/ && rm ./$input_folder/GRCh38.d1.vd1.fa.tar.gz
#gtf
echo "PREPARE gtf file\n"
wget -O ./$input_folder/gencode.v36.annotation.gtf.gz https://api.gdc.cancer.gov/data/be002a2c-3b27-43f3-9e0f-fd47db92a6b5
gunzip -c ./$input_folder/gencode.v36.annotation.gtf.gz > ./$input_folder/gencode.v36.annotation.gtf && rm ./$input_folder/gencode.v36.annotation.gtf.gz
#pon
echo "PREPARE PON file\n"
wget -O ./$input_folder/1000g_pon.hg38.vcf.gz https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
wget -O ./$input_folder/1000g_pon.hg38.vcf.gz.tbi https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi
#germline
echo "PREPARE germline resource file\n"
wget -O ./$input_folder/af-only-gnomad.hg38.vcf.gz https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
wget -O ./$input_folder/af-only-gnomad.hg38.vcf.gz.tbi https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi
#RNA-EDI
echo "PREPARE RNA-edits resource file\n"
wget -O ./$input_folder/RNA-EDI.txt.gz http://srv00.recas.ba.infn.it/webshare/ATLAS/donwload/TABLE1_hg38.txt.gz
gunzip -c ./$input_folder/RNA-EDI.txt.gz > ./$input_folder/RNA-EDI.txt && rm ./$input_folder/RNA-EDI.txt.gz
awk 'NR>1 {print $1":"$2"-"$2}' ./$input_folder/RNA-EDI.txt > ./$input_folder/RNA-EDI_38.bed && rm ./$input_folder/RNA-EDI.txt
#dbsnp
echo "PREPARE dbSNP file\n"
wget -O ./$input_folder/common_all_20180418.vcf.gz https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/GATK/common_all_20180418.vcf.gz
wget -O ./$input_folder/common_all_20180418.vcf.gz.tbi https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/GATK/common_all_20180418.vcf.gz.tbi
#gene_list
echo "PREPARE gene list file\n"
perl ./tools/merge_genes.pl ./$input_folder/gencode.v36.annotation.gtf ./$input_folder/gg.list
cat ./$input_folder/gencode.v36.annotation.gtf | grep $'\tgene\t' | grep 'HLA-' > ./$input_folder/HLA.gtf
cat ./$input_folder/gencode.v36.annotation.gtf | grep $'\tgene\t' | grep 'IG_' > ./$input_folder/IGG.gtf
cat ./$input_folder/gencode.v36.annotation.gtf | grep $'\tgene\t' | grep 'pseudogene' > ./$input_folder/pseudoGene.gtf

#generate indices
./$tools_folder/samtools-1.10/samtools faidx ./$input_folder/GRCh38.d1.vd1.fa
java -jar ./$tools_folder/picard.jar CreateSequenceDictionary -R ./$input_folder/GRCh38.d1.vd1.fa -O ./$input_folder/GRCh38.d1.vd1.dict
./$tools_folder/hisat2-2.1.0/hisat2-build -p 8 ./$input_folder/GRCh38.d1.vd1.fa ./$input_folder/hisat2_GRCh38