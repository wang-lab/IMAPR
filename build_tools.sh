#!/bin/bash

input_folder="$1"

if [ -z "$input_folder" ]; then
  echo "Error: No output tools folder were provided, please specifiy where you want to store the tools file."
  exit 1
fi

if [ ! -d "$input_folder" ]; then
	mkdir $input_folder
fi

echo "Download hisat2\n"
if [ -d "./$input_folder/hisat2-2.1.0" ]; then
	rm -r ./$input_folder/hisat2-2.1.0
fi
wget -O ./$input_folder/hisat2-2.1.0-Linux_x86_64.zip https://cloud.biohpc.swmed.edu/index.php/s/hisat2-210-Linux_x86_64/download
unzip ./$input_folder/hisat2-2.1.0-Linux_x86_64.zip -d ./$input_folder && rm ./$input_folder/hisat2-2.1.0-Linux_x86_64.zip

echo "Download gatk tool kit\n"
if [ -d "./$input_folder/gatk-4.1.8.1" ]; then
	rm -r ./$input_folder/gatk-4.1.8.1
fi
wget -O ./$input_folder/gatk-4.1.8.1.zip https://github.com/broadinstitute/gatk/releases/download/4.1.8.1/gatk-4.1.8.1.zip
unzip ./$input_folder/gatk-4.1.8.1.zip -d ./$input_folder && rm ./$input_folder/gatk-4.1.8.1.zip
wget -O ./$input_folder/picard.jar https://github.com/broadinstitute/picard/releases/download/2.23.5/picard.jar

echo "Download samtools kit\n"
if [ -d "./$input_folder/samtools-1.10" ]; then
	rm -r ./$input_folder/samtools-1.10
fi
wget -O ./$input_folder/samtools-1.10.tar.bz2 https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
tar -xvjf ./$input_folder/samtools-1.10.tar.bz2 -C ./$input_folder && rm ./$input_folder/samtools-1.10.tar.bz2
cd ./$input_folder/samtools-1.10
autoheader
autoconf -Wno-syntax
./configure
make

cd ..
cd ..
if [ -d "./$input_folder/bcftools-1.10.2" ]; then
	rm -r ./$input_folder/bcftools-1.10.2
fi
wget -O ./$input_folder/bcftools-1.10.2.tar.bz2 https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2
tar -xvjf ./$input_folder/bcftools-1.10.2.tar.bz2 -C ./$input_folder && rm ./$input_folder/bcftools-1.10.2.tar.bz2
cd ./$input_folder/bcftools-1.10.2
./configure
make