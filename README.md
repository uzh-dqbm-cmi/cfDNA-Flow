# cfDNA-Flow


## 1. Overview
cfDNA-Flow facilitates the accurate and reproducible analysis of cfDNA WGS data. It offers various preprocessing options to accommodate different experimental setups and research needs in the field of liquid biopsies. 

![](/Users/ivna/00_projects/cfDNA-Flow/workflow.png)

## 2. Preprocessing options
### 2.1 Trimming Options
cfDNA-Flow provides the flexibility to either trim or not trim the input reads based on the user's requirements. Trimming removes low-quality bases, which can impact downstream analyses.

### 2.2 Reference Genome Selection
Users can choose from the following genome builds: hg19, hg38, hg19decoy, and hg38noalt.

### 2.3 Post-Alignment Filtering and GC bias correction
The pipeline utilizes the BWA-MEM software for alignment, followed by extensive post-alignment filtering steps to ensure reliable alignments. Users can define specific filtering criteria to remove low-quality or ambiguous reads, such as secondary alignments, reads with insertion or deletion, and reads with low mapping qualities. Additionally, the pipeline offers an option to correct for GC bias.

## 3. Usage
To use cfDNA-Flow, follow these steps:

### 3.1 Installation:
        Clone the cfDNA-Flow repository from GitHub.
        Install the required dependencies...
### 3.2 Configuration:
        Edit the configuration file to specify the input files, reference genome, and desired preprocessing options.
        Customize the trimming, alignment, and filtering settings as needed.
### 3.3 Execution/Demo:
        Run the pipeline...
### 3.4 Output:
        The pipeline outputs processed reads, alignment files (BAM files, BED files), and comprehensive quality control reports.

## 4. Support
For issues, questions, or contributions, please visit the cfDNA-Flow GitHub repository or contact the maintainers via email. Detailed documentation, including installation guides, is available in the repository's wiki.
