# IMAPR
This repository holds the dual alignment based RNA-Seq somatic mutation, IMAPR, pipeline used to discover RNA-seq based somatic mutation.

![workflow](https://github.com/wang-lab/IMAPR/blob/main/workflow.PNG)

The pipeline accepts aligned RNA-seq data as inputs and identifies somatic mutations within it.
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
   * Note: higher versions of GATK may not be compatible with IMAPR due to potential deprecation or changes to parameters used in GATK.
   
* The Picard tool v2.22.8 is required.
   * [Picard download](https://github.com/broadinstitute/picard/releases)
   * Note: higher versions of Picard may not be compatible with IMAPR due to potential deprecation or changes to parameters used in Picard.
   
* The samtools version 1.10 or higher is required.
   * [samtools download](https://github.com/samtools/samtools)
   
   
### Prerequisites

Follwing reference files are required to run IMAPR pipeline:

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
	* Three RNA-edits resource files were provided in IMAPR/reference
	* IMAPR/reference/Darned_38.bed
	* IMAPR/reference/Radar_38.bed
	* IMAPR/reference/DRNA-EDI.bed
	
### Installation of IMAPR standalone program

* Install IMAPR using git command:
```
   git clone https://github.com/wang-lab/IMAPR.git   
```
* Prepare your input files and update your input.txt file.
* Type 'sh IMAPR.sh' to run the program and view the help file.

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
  hisat2	<path-to-hisat2-package>

  #reference
  fasta_ref	<path-to-fasta-reference>
  gtf_ref	<path-to-gtf-reference>
  genelist_ref	<path-to-gene-list-reference>  
  dbsnp_ref	<path-to-dbsnp-reference>
  germline_ref	<path-to-germline-reference>
  PON_ref	<path-to-PON-reference>
  hisat_ref	<path-to-hisat-ref-reference>
  igg_ref	<path-to-igg-reference>
  hla_ref	<path-to-hla-reference>
  pseudo_ref	<path-to-pseudo-reference>
  tcga_PON_ref	<path-to-tcga-PON-reference>
  radar_ref	<path-to-radar-reference>
  darned_ref	<path-to-darned-reference>
  REDI_ref	<path-to-REDI-reference>  
```  
***sample_name*** is the name of the sample, used to name the output folder of the run. For consistency we recommend that sample_name is the same as the name of the raw data folder.  
***input_format*** dictates whether the input files both RNA-seq data or RNA vs WXS data.  
***tumor_input*** is path to aligned bam file for tumor sample.  
***normal_input*** is path to aligned bam file for normal sample.  
***out_prefix*** is path to output folder.  

***gatk*** is path to gatk tool package  
***picard*** is path to picard jar file.  
***samtools*** is path to samtools package.  
***hisat2*** is path to hisat2 package.  

***fasta_ref*** is path to genome fasta reference.  
***gtf_ref*** is path to gtf reference.  
***genelist_ref*** is path to gene list reference.  
***RNA_edits_ref*** is path to RNA-edits reference.  
***dbsnp_ref*** is path to dbsnp reference.  
***germline_ref*** is path to germline reference.  
***PON_ref*** is path to PON reference.  
***hisat_ref*** is path to hisat reference.  
***igg_ref*** is path to igg gene list reference.  
***hla_ref*** is path to hla gene list reference.  
***pseudo_ref*** is path to pseudo-gene list reference.  
***tcga_PON_ref*** is path to tcga PON reference.  
***radar_ref*** is path to RADAR reference.  
***darned_ref*** is path to DARNED reference.  
***REDI_ref*** is path to REDI reference.  

### Command Line Parameters

* bash command submission, user need to change the indices listed in IMAPR.sh   
   
```
   sh IMAPR.sh pipeline_inputs.txt
```

* The script can also be run separately with follow steps
  
```
   perl detect_variants.pl [options]... -ID $sample_name -mode $input_format -T $tumor_input -N $normal_input -R $fasta_ref -O $out_prefix -gatk $gatk -picard $picard -hisat2 $hisat2 -gtf $gtf_ref -gene $genelist_ref -dbsnp $dbsnp_ref -hisat2_reference $hisat_ref -germline $germline_ref -pon $PON_ref
   perl filter_variants.pl [options]... -ID $sample_name -O $out_prefix -R $fasta_ref -igg $igg_ref -hla $hla_ref -pseudo $pseudo_ref -tcga $tcga_PON_ref -radar $radar_ref -darned $darned_ref -redi $REDI_ref
   perl machine_learning.pl [options]... -ID $sample_name -O $out_prefix -gtf $gtf_ref
```

### Outputs

If the file is read in correctly, the following output files will be generated in output folder.
* 1st Mutect2 Variants files
```
   [sample_ID]_first_variants.txt
```
* Repeated Mutect2 Variants files
```
   [sample_ID]_final_variants.txt
```
* samtools mpileup Variants files
```
   [sample_ID]_bcftools.output
```
* Machine-learning filtered Variants files
```
   [sample_ID]_mc_filter_Variants.vcf
```
## License & copyright

License under the [GNU public library](LICENSE)