---
output:
  pdf_document: default
  html_document: default
---
# Secondary Findings Tool

## Overview

This repository contains a novel tool designed for the automated management of secondary findings. It enables users to analyze VCF (Variant Call Format) files for automatic management of secondary findings related to personal risk, reproductive risk, and pharmacogenetics.

## <a name="TOC">Table of content</a>
 * [Installation](#installation)
 * [Usage](#usage)
 * [Outputs](#outputs)
 * [Example](#example)
 * [Panel creation](#panel) 
 * [Data](#data)  
 * [Dependencies](#dependencies)
 * [Citation](#citation) 
 * [Version history](#versionhistory)


## <a name="installation">Installation</a>

Before using this tool, please ensure you have the following prerequisites:
1. Linux Environment: This tool is designed to work in a Linux environment. If you're using Windows, consider using Windows Subsystem for Linux (WSL) or a Linux virtual machine.
2. Python: Make sure you have Python 3.8 or higher installed. You can check your Python version by running:

```
python --version
```
If Python is not installed, you can download it from the official Python website (https://www.python.org/downloads/) or use your system's package manager.

3. Git: You'll need Git to clone the repository. If Git is not already installed, you can install it using your system's package manager or download it from the official Git website (https://git-scm.com/downloads).

4. InterVar and ANNOVAR: These are external tools that this software relies on. Please ensure that you have both InterVar (https://github.com/WGLab/InterVar) and ANNOVAR (https://annovar.openbioinformatics.org/en/latest/) installed in your system. You can find installation instructions for these tools on their respective websites.

5. Also, this tool uses hs37d5 as reference sequence, which is available at:  https://console.cloud.google.com/storage/browser/genomics-public-data/references/hs37d5?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false.

Once you have the prerequisites in place, you can proceed with the installation:
1. Clone this repository to your local machine using Git:
```
git clone https://github.com/edurlaf/tfm_secondaryfindings.git
```
2. Navigate to the project directory:
```
cd tfm_secondaryfindings
```
3. Now, you need to configure your tool by providing the necessary parameters. These parameters are stored in a configuration file named "config.json". You should customize this file according to your specific environment and requirements.
4. After configuring the tool, you should be ready to use it as described in the "Usage" section of this README.
## <a name="Usage">Usage</a>

This tool is designed to assist in the automated management of secondary findings in genomic data. You can run it from the command line with various options to customize its behavior. The basic command structure is as follows:


```
python secondary_findings.py input_file.vcf \
         [--mode <Option: 'basic' or 'advanced'>]
         [--evidence <integer>]
         [--assembly <Option: '37' or '38'>]
```
where:
 * **input_file.vcf**: This should be replaced with the path to your VCF (Variant Call Format) file that contains the genomic data you want to analyze. 
 * **--mode**: (Optional) This option allows you to choose the analysis mode. You can select either 'basic' (default) or 'advanced' based on your requirements. The 'basic' mode runs InterVar, while the 'advanced' mode combines InterVar with ClinVar database.
 * **--evidence**: (Optional) Use this option to specify the evidence level or review status you want to apply to ClinVar database. Provide an integer value from 1-4 to set the evidence level. Default: 1.
 * **--assembly**: (Optional) Select the genome assembly version you want to use for the analysis. You can choose either '37' (default) or '38' depending on the assembly that corresponds to your data.



## <a name="outputs">Outputs</a>

After running the tool, you will find various output files that summarize the analysis of secondary findings in genomic data. These outputs are generated in the designated folders.
1. **Final output**: The main result is stored in an Excel file, located in the "final_output" directory. This Excel file contains several sheets:
* **Personal Risk**: this sheet presents secondary findings related to personal risk.
* **Reproductive Risk**: this sheet compiles findings associated with reproductive risk.
* **Pharmacogenetic Risk**: here, you can find the results for pharmacogenetic risk.
* **Diplotypes and Phenotypes**: this fourth sheet contains information on diplotypes, phenotypes, and activity scores for 5 pharmacogenetic risk genes.

2. **Intermediate Outputs** (in the "temp" directory):
    * normalized.vcf: This file stores the results after the normalization process of the input VCF file.
    * intersection.vcf: It contains the findings after the intersection of VCF data with predefined BED files for each category.
    * multianno.intervar: You can locate the outcomes of the InterVar tool's analysis here.
    * all_results.csv: In this CSV file, you will find all pathogenic (P) or likely pathogenic (LP) variants before filtering based on inheritance rules.

## <a name="output">Example</a>

An easy example is provided for testing purposes. To test this example you can just run (after [Installation](#installation)):

```
# Docker mode
./impuSARS --infile example/sequence.fa \
           --outprefix imputation 
# Conda mode
conda activate impusars
impuSARS --infile example/sequence.fa \
         --outprefix imputation
conda deactivate 
```

The [example SARS-CoV-2 sequence](example/sequence.fa) has been internally sequenced and is available under the ENA Accession [PRJEB43882](https://www.ebi.ac.uk/ena/browser/view/PRJEB43882) (see [Data](#data) for details). This sequence includes a high rate of missing regions (Ns). Therefore, impuSARS will return a completely imputed genome sequence (FASTA file) and its corresponding assigned lineage (CSV file).

## <a name="panel">Panel creation</a>

impuSARS tool now includes another all-in script for users to create their own reference panel for SARS-CoV-2 or any other viral sequences to impute. Reference panels can be created as follows:

```
# Docker mode
./impuSARS_reference --name <reference_prefix> \
                     --output_path <output_path> \
                     --input_fasta <input_fasta> \
                     --genome_fasta <reference_fasta> \
                     [--unknown_nn <unknown_nn>]
                     [--threads <num_threads>] 
# Conda mode
conda activate impusars
impuSARS_reference --name <reference_prefix> \
                     --output_path <output_path> \
                     --input_fasta <input_fasta> \
                     --genome_fasta <reference_fasta> \
                     [--unknown_nn <unknown_nn>]
                     [--threads <num_threads>]
conda deactivate
```
where:
 * **<output_path>**: Directory where the custom reference panel will be generated.
 * **<reference_prefix>**: prefix name given to the output reference panel without extension. Output will generate <reference_prefix>.m3vcf.gz reference panel file.
 * **<input_fasta>**: FASTA file including the alignment of all sequences used to train and generate the reference panel.
 * **<genome_fasta>**: FASTA file with the reference genome for the virus to impute. For example, [SARS-CoV-2 reference](docker_files/references/SARS_CoV_2_REFERENCE.v1.0.fasta).
 * **<unknown_nn>**: (Optional) Special character used in alignment for missing nucleotides, if any. Default: "n".
 * **<num_threads>**: (Optional) Number of CPUs used for imputation. Default: 1.

As before, experienced users can run the script directly using Docker as:

```
docker run -it --rm -v <input_path>:/data -v <ref_path>:/ref -v <output_path>:/output impusars \
       impuSARS_reference --name <reference_prefix> \
                          --output_path /output/ \
                          --input_fasta /data/<input_fasta_basename> \
                          --genome_fasta /ref/<genome_fasta_basename> \
                          [--unknown_nn ${unknn}] \
                          [--threads ${threads}]
```
where **<input_path>, <ref_path>** refer to directories where <input_fasta> and <genome_fasta> are respectively located whereas **<input_fasta_basename>** and **<genome_fasta_basename>** are the basenames of those files (without path). 

## <a name="output">Data</a>

Nine internally sequenced SARS-CoV-2 samples are available at the following repository for validation purposes:

* **Raw sequencing data and consensus sequences:**: [ENA Dataset Accession ID PRJEB43882](https://www.ebi.ac.uk/ena/browser/view/PRJEB43882).
* **ImpuSARS imputed sequences and lineages:**: [Zenodo repository](https://doi.org/10.5281/zenodo.4616731).

Also, impuSARS uses the [hCoV-19/Wuhan/WIV04/2019](https://www.ncbi.nlm.nih.gov/nuccore/MN908947) sequence as the official reference sequence, which is available [here](docker_files/references/SARS_CoV_2_REFERENCE.v1.0.fasta).

Finally, impuSARS was initially trained with a reference panel containing 239,301 sequences from [GISAID](https://www.gisaid.org/) (downloaded by January 7, 2021). Therefore, we would like to gratefully acknowledge all those laboratories and sequence contributors that made possible to create such a reference panel ([acknowledgment](acknowledgement/gisaid_hcov-19_acknowledgement_table_2021_04_27_10.pdf)). **Current reference version (v2.1) contains 899,447 sequences (updated by June 17th, 2021).**

## <a name="dependencies">Dependencies</a>

This tool relies on the following dependencies:

 * [InterVar](https://github.com/samtools/bcftools) (v1.11)
 * [Annovar](https://www.drive5.com/muscle/) (v3.8.31)

Since impuSARS is encapsulated in a Docker image to facilitate distribution, **only Docker installation is required**. Docker can be downloaded for any operating system at [Get Docker](https://docs.docker.com/get-docker/). In case **conda installation** is preferred, please note that two command packages are required:

 * [Conda](https://docs.conda.io/en/latest/) 
 * **curl** or **wget** for downloading dependencies.


Date

    August 1, 2023

Contact

    Email: edurlaf@gmail.com
    GitHub: github.com/edurlaf