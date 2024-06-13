# cfDNA-Flow


## 1. Overview
cfDNA-Flow facilitates the accurate and reproducible analysis of cfDNA WGS data. It offers various preprocessing options to accommodate different experimental setups and research needs in the field of liquid biopsies. 

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

### 3.3 Normalized coverage
cfDNA-Flow splits the genome into 1-Mbp bins and counts the number of fragments in each region. It then scales the 1-Mbp bin-wise fragment counts by dividing them by the average number of fragments across all bins. This scaled coverage is used to calculate Pearson’s correlation between healthy and cancer samples. When calculating the correlation of the coverage of a healthy sample to the average coverage of healthy samples, the given healthy sample is excluded from the average.

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
        Edit the configuration file to specify the input files, reference genome, and desired preprocessing options.
        Customize the trimming, alignment, and filtering settings as needed.

### 4.3 Execution/Demo:
To start preprocessing, execute the following command. Use the -np flag for a dry run to verify everything works correctly.

        snakemake -s Snakefile --configfile test/test_cfDNA_pipeline.yaml -j 2 do_preprocess -np 

If successful, rerun the command without the -np flag.

        snakemake -s Snakefile --configfile test/test_cfDNA_pipeline.yaml -j 2 do_preprocess 

Next, bed to process BED files, use:

        snakemake -s Snakefile --configfile test/test_cfDNA_pipeline.yaml -j 2 do_bedprocess
        
Do unique counts:

        snakemake -s Snakefile --configfile test/test_cfDNA_pipeline.yaml -j 2 do_uniq_counts

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
The pipeline outputs processed reads, alignment files (BAM files, BED files), and comprehensive quality control reports.

## 5. Support
For issues, questions, or contributions, please visit the cfDNA-Flow GitHub repository or contact the maintainers via email. 