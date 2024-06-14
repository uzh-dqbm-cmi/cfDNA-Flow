# cfDNA-Flow


## 1. Overview
cfDNA-Flow facilitates the accurate and reproducible analysis of cfDNA WGS data. It offers various preprocessing options to accommodate different experimental setups and research needs in the field of liquid biopsies. 

![](https://github.com/uzh-dqbm-cmi/cfDNA-Flow/blob/main/workflow.png)

## 2. Preprocessing options
### 2.1 Trimming Options
cfDNA-Flow provides the flexibility to either trim or not trim the input reads based on the user's requirements. Trimming removes low-quality bases, which can impact downstream analyses.

### 2.2 Reference Genome Selection
Users can choose from the following genome builds: hg19, hs37d5 (hg19decoy), hg38 and hg38 without alternative contigs (hg38noalt). For download links, please refer to the Data Availability section in the accompanying paper.

### 2.3 Post-Alignment Filtering and GC bias correction
The pipeline uses the BWA software for alignment, followed by extensive post-alignment filtering steps to ensure reliable alignments. Users can define specific filtering criteria to remove low-quality or ambiguous reads, such as secondary alignments, reads with insertion or deletion, and reads with low mapping qualities. Additionally, the pipeline offers an option to correct for GC bias.

## 3. Feature Extraction

### 3.1 Fragment length features
cfDNA-Flow offers fragment length analysis; calculating the mean, median, and standard deviation values for fragments sized 100 to 220 base pairs (bp), corresponding to the mononucleosomal size range. Additionally, cfDNA-Flow calculates the frequencies of cfDNA fragment sizes ranging from 70 bp to 1000 bp in 10 bp bins.

### 3.2 Copy number changes
cfDNA-Flow utilizes two copy number analysis tools: ichorCNA (v0.2.0) and tMAD, to estimate copy number changes and tumor fraction.

### 3.3 Binned coverage 
cfDNA-Flow splits the genome into 1-Mbp bins and counts the number of cfDNA fragments in each region. 

### 3.4 Differential coverage analysis over DNase hypersensitivity sites
This feature was calculated using [LIQUORICE](https://github.com/epigen/LIQUORICE/tree/master). Please note that the cfDNA-Flow pipeline does not integrate the LIQUORICE tool for detecting epigenetic signatures in cell-free DNA from liquid biopsies. However, you can find more information and access LIQUORICE through the following [link](https://liquorice.readthedocs.io/en/latest/).

## 4. Usage
To use the cfDNA-Flow, follow these steps:

### 4.1 Installation:
Clone the cfDNA-Flow repository from GitHub to your local machine.

        git clone https://github.com/uzh-dqbm-cmi/cfDNA-Flow.git
        cd cfDNA-Flow

It is recommended to create a virtual environment to manage the project's dependencies. This ensures that the dependencies do not interfere with other Python projects on your machine. See how to create a Python environment [here](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/).

Once the virtual environment is activated, install the required Python dependencies using the `requirements.txt` file.

        pip install -r requirements.txt

Additionally, some R packages are required for the project. Make sure you have R installed (version 4.3). You can install these packages by running the script below:

        R -f install_packages.R

After following these steps, your environment should be set up with all the necessary dependencies for both Python and R. You are now ready to proceed with using the cfDNA-Flow pipeline. See section 4. Usage. 

Once you are finished using cfDNA-Flow, deactivate the virtual environment by running:

        deactivate

### 4.2 Configuration:
The configuration file, `test_cfDNA_pipeline.yaml`, is used to specify the input files, reference genome, and desired preprocessing options. User can customize the trimming, alignment, and filtering settings as needed.
Settings for this demo are as follows: reads are trimmed, reference genome is hg38, mapping quality is 30, [SAM flag](https://broadinstitute.github.io/picard/explain-flags.html) is 256, CIGAR string is D. 

### 4.3 Execution/Demo:
To start preprocessing, execute the following command. Use the -np flag for a dry run to verify everything works correctly.

        snakemake -s Snakefile --configfile test/test_cfDNA_pipeline.yaml -j 2 do_preprocess -np 

If successful, rerun the command without the -np flag.

        snakemake -s Snakefile --configfile test/test_cfDNA_pipeline.yaml -j 2 do_preprocess 

Next, bed to process BED files, use:

        snakemake -s Snakefile --configfile test/test_cfDNA_pipeline.yaml -j 2 do_bedprocess

Do global length: 

        snakemake -s Snakefile --configfile test/test_cfDNA_pipeline.yaml -j 2 do_global_length

Do tMAD:

        snakemake -s Snakefile --configfile test/test_cfDNA_pipeline.yaml -j 2 do_cal_blacklist
        snakemake -s Snakefile --configfile test/test_cfDNA_pipeline.yaml -j 2 do_cal_RefSample
        snakemake -s Snakefile --configfile test/test_cfDNA_pipeline.yaml -j 2 do_cal_t_MAD_forall
        snakemake -s Snakefile --configfile test/test_cfDNA_pipeline.yaml -j 2 do_visualising_t_MAD_forall

Do ichorCNA:

        snakemake -s Snakefile --configfile test/test_cfDNA_pipeline.yaml -j 2 do_createPoN
        snakemake -s Snakefile --configfile test/test_cfDNA_pipeline.yaml -j 2 do_ichorCNA
        snakemake -s Snakefile --configfile test/test_cfDNA_pipeline.yaml -j 2 do_ichorCNA_results

### 4.4 Output:
The pipeline outputs alignment files (BAM files, BED files), and quality control reports. Additionally, it outputs features described in section 3. Feature Extraction. 

#### BAM files
Processed BAM files of studied samples, accompanied by their .bai index files, are stored in the `results/BAM/FalseD25630` folder and have `.sortByCoord.bam` suffix. 

#### BED files
BED files are stored in the `results/BED/FalseD25630` folder. Those files store information about chromosome number, start and end positions of cfDNA fragments.

#### Quality control reports
Output of QC is stored in the `results/QC/FalseD25630/multiqc_data` folder. Specifically, `multiqc_report.html` file contains multiple the QC metrics: general statistics, Picard metrics (alignment summary, mean read length, mark duplicates, WGS coverage, WGS filtered bases), FastQC (sequence counts, sequence quality histograms and quality scores, per base sequence content, per sequence GC content, per base N content, sequence length distribution, sequence duplication levels, overrepresented sequences, adapter content, status checks).  

#### Fragment length features
The output of fragment length features is stored in the `results/feature/FalseD25630/global_length.tsv` file. Columns store fragment length features for each studied sample (rows).

#### Coverage features and fragment lengths in 1Mb genomic bins
Outputs of features in 1Mb genomic bins can be found in the `results/BED/FalseD25630` folder. Values for all the samples are stored in `medgeddf.csv` file. Values for each individual sample are stored in the files with suffix `binned.csv`.

Additional length features for every sample are stored in the folder `results/BED/FalseD25630` and have the following suffixes:

`binned_lengths.csv` - each row contains information about the chromosome number, genomic bin number (1Mb wide), and the lengths of all cfDNA fragments corresponding to that bin

`len.csv` - contains a single column listing the lengths of all cfDNA fragments derived from a sample

`lenuniqcount.csv` - a two-column format representing the histogram of cfDNA fragment lengths along with their frequencies

#### ichorCNA
Results of ichorCNA analysis can be found in the `results/feature/FalseD25630/ichorCNA` folder. For detailed ichorCNA output description see this [link](https://github.com/broadinstitute/ichorCNA/wiki/Output). Shortly, ichorCNA outputs tumor fraction estimates based on CNA analysis. Additionally, it outputs CNA plots representing log2 ratio copy number for each bin in the genome.

![](https://github.com/uzh-dqbm-cmi/cfDNA-Flow/blob/main/ichorCNA_plot.png)

#### tMAD
The outputs of [tMAD analysis](https://github.com/sdchandra/tMAD) are stored in the `results/BED/FalseD25630/tMAD/tMAD_results.tsv` file. This file contains sample names and their corresponding tMAD values. 

## 5. Support
For issues, questions, or contributions, please visit the cfDNA-Flow GitHub repository or contact the maintainers via email. 