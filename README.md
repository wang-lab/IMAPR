# DARSM
This repository holds the dual alignment based RNA-Seq somatic mutation, DARSM, pipeline used to discover RNA-seq based somatic mutation.


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

* This package is supported for *Linux* operating systems.  The package has been tested on the following systems:
```
   Linux: Ubuntu 18.04.4 LTS
```
* Perl 5 interpreter or higher on a Ubuntu compatible Linux system is required.
   * [Installation instruction](https://learn.perl.org/installing/)
   * [Perl download](https://www.perl.org/get.html)
   
* Python 3.8.10 or higher on a Ubuntu compatible Linux system is required.   
   * [Python3 download](https://www.python.org/downloads/)
   
* Following python3 packages are required.
   * sklearn
   * joblib
   * pandas
   * numpy
   
* The Genome Analysis Toolkit (GATK) v4.1.8.1 is required.
   * [GATK download](https://github.com/broadinstitute/gatk/releases)
   * Note: higher versions of GATK may not be compatible with DARSM due to potential deprecation or changes to parameters used in GATK.
   
* The Picard tool v2.22.8 is required.
   * [Picard download](https://github.com/broadinstitute/picard/releases)
   * Note: higher versions of Picard may not be compatible with DARSM due to potential deprecation or changes to parameters used in Picard.
   
* The samtools version 1.10 or higher is required.
   * [samtools download](https://github.com/samtools/samtools)
   
   
### Prerequisites

Follwing reference files are required to run DARSM pipeline:

* A fasta reference sequence file:
	* The fasta reference sequence, GRCh38.d1.vd1.fa.tar.gz, from TCGA data portal is recommended(https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files)
	
* A annotation files reference sequence file:
	* The annotation files from TCGA data portal is recommended(https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files)
	
* A panel of normal (PON) file:
	* The annotation files from TCGA data portal is recommended(https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files)
	* The PON files recommeded by broad institute is recommended(https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38)
	
* A Germline resource file:
	* The annotation files from TCGA data portal is recommended(https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38)

* Three RNA-edits resource files:
	* Three RNA-edits resource files were provided in DARSM/reference
	* DARSM/reference/Darned_38.bed
	* DARSM/reference/Radar_38.bed
	* DARSM/reference/DRNA-EDI.bed
	
### Installation of DARSM standalone program

* Install DARSM using git command:
```
   git clone https://github.com/wang-lab/DARSM.git   
```
* Prepare your input files and update your input.txt file.
* Type 'sh DARSM.sh' to run the program and view the help file.

### I/O Descriptions  
#### Inputs  
***pipeline_inputs.txt***  
This file has the following format. The order of rows doesn't matter.  

```
  sample_name	<sample name>
  input_format	<RNA/RNA OR RNA/WXS>
  tumor_input	<path-to-bam-file>
  normal_input	<path-to-bam-file>
  out_prefix	<path-to-output-folder>
  
  #tools_reference
  gatk	<path-to-gatk-package>
  picard	<path-to-picard-package>
  samtools	<path-to-samtools-package>
  
  #reference
  fasta_ref	<path-to-fasta-reference>
  RNA_edits_ref	<path-to-rna-edits-referenc>
  germline_ref	<path-to-germline-referenc>
  PON_ref	<path-to-PON-referenc>
```  
***sample_name*** is the name of the sample, used to name the output folder of the run. For consistency we recommend that sample_name is the same as the name of the raw data folder  
***input_format*** dictates whether the input files both RNA-seq data or RNA vs WXS data.  
***tumor_input*** is path to aligned bam file for tumor sample.  
***normal_input*** is path to aligned bam file for normal sample. 
***out_prefix*** is path to output folder. 

***gatk*** dictates whether the input files both RNA-seq data or RNA vs WXS data.  
***picard*** is path to aligned bam file for tumor sample.  
***samtools*** is path to aligned bam file for normal sample. 

***fasta_ref*** is path to genome fasta reference.  
***RNA_edits_ref*** is path to RNA-edits reference.  
***germline_ref*** is path to germline resource reference.  
***PON_ref*** is path to panel of normal reference.  

### Command Line Parameters

* bash command submission, user need to change the indices listed in DARSM.sh   
   
```
   sh DARSM.sh pipeline_inputs.txt
```

* The script can also be run separately with follow steps
  
```
   perl p1.pl [options]... -Tumor *.bam -Normal *.bam -R /directory/to/reference -outPrefix /directory/to/output
   perl p2.pl [options]... -inputDir *.bam -/directory/to/p1.pl/output
   perl p3.pl [options]... -inputDir *.bam -/directory/to/p1.pl/output
```

### Outputs

If the file is read in correctly, the following output files will be generated in output folder.
* 1st Variants files
```
   first_variants.txt
```
* Variants files
```
   final_variants.txt
```
* Machine-learning filtered Variants files
```
   machine_learning_variants.txt
```
* Filtered somatic mutation files
```
   filtered_final_variants.txt
```
* annotated somatic mutation files
```
   anno_filtered_final_variants.txt
```
## License & copyright

License under the [GNU public library](LICENSE)