# iPrimer
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

### Installation of iPrimer standalone program

* Place the DARSM.tar.gz file anywhere in your Linux system and uncompress using the following command:
```
   tar -xzvf DARSM.tar.gz   
```
* Copy your input files into the newly created DARSM directory.
* Type 'sh DARSM.sh' to run the program and view the help file.

### Command Line Parameters

* Bam file submission, required (-f path/to/fasta).
   This option allows the user to submit one sequence in a FASTA file, using the following command:
```
   sh DARSM.sh
```

### Outputs

If the file is read in correctly, the following output files will be generated in iPrimer folder.

* Variants files
```
   final_variants.txt
```
* Filtered somatic mutation files
```
   filtered_final_variants.txt
```
* Filtered mutation files
```
   anno_filtered_final_variants.txt
```
## License & copyright

License under the [GNU public library](LICENSE)