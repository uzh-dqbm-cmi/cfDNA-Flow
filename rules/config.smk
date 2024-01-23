import sys
from distutils.util import strtobool

global WORKDIR
global GENOMEBIT
global REFGENOME
global PARAMDIR
global SAMPLEDIR
global ANNO_NAME
global ANNO_FILE
global shendurelab_cfDNA
global ANNO_FILE
global SIZE_SELECTION_MIN
global SIZE_SELECTION_MAX
global SIZE_SELECTION_SUB_FOLDER
global BIN_SIZE
global CHROMS
global SCRIPTS
global THREADS
global TEMP
global SEED

# still hard coded in some places
global GENOME
global CNA_DO_CONTROL_REFS
global RC_NORMALIZATION
global hash_it
global MIXED_MODEL_FILE_SUFF
global WPS_WINDOW
global WPS_MIN
global WPS_MAX
global WPS_BIGWIG_NORMALIZED
global WPS_BIN_SIZE
global WPS_BIN_NORMALIZATION
global KEEP_BEDGRAPH
global FFT_BINSIZE
global FFT_FREQ_FROM
global FFT_FREQ_TO
global FFT_FREQ_SEL


WORKDIR = config["outputDir"]
GENOMEBIT = config["genomeBit"]
REFGENOME = config["refGenome"]
PARAMDIR = config["paramDir"]
SAMPLEDIR = config["sampleDir"]
ANNO_NAME = config["anno_name"]
ANNO_FILE = config["anno_file"]
shendurelab_cfDNA = config["shendurelab_cfDNA"]
cfNOMe = config["cfNOMe"]
SIZE_SELECTION_MIN = str(config["size_selection"]["min"]) if config["size_selection"]["min"] else ""
SIZE_SELECTION_MAX = str(config["size_selection"]["max"]) if config["size_selection"]["max"] else ""
# scripts = config["scripts"]
THREADS = config["threads"]
SCRIPTS = config["scripts"]
CONTROL_SAMPLES = config["control_samples"]   # file name: every line is sample name
SAMPLES = config["samples"]   # file name: every line is sample name
try:
    TEMP = config["temp"]   # temporary folder
except:
#    print("Missing 'temp' configuration parameter. error:", sys.exc_info()[0])
    TEMP = ""   # used defaults

try:
    SEED = config["seed"]   # Set the seed of the random number generator, if you like to produce the same results.
except:
#    print("Missing 'seed' configuration parameter. error:", sys.exc_info()[0])
    SEED = ""   # == don't set seed

BIN_SIZE = str(config["copy_number_aberration"]["bin_size"])   # "1000"
WPS_WINDOW = str(config["wps"]["window"]) if config["wps"]["window"] else "120" # default according https://doi.org/10.1016/j.cell.2015.11.050
WPS_MIN = str(config["wps"]["min"]) if config["wps"]["min"] else "120"  # default min fragments according https://doi.org/10.1016/j.cell.2015.11.050
WPS_MAX = str(config["wps"]["max"]) if config["wps"]["max"] else "180"  # default max fragments according https://doi.org/10.1016/j.cell.2015.11.050
WPS_BIN_SIZE = str(config["wps"]["bin_size"]) if config["wps"]["bin_size"] else "1000"
try:
    WPS_BIN_NORMALIZATION = str(config["wps"]["bin_normalization"])
except:
    WPS_BIN_NORMALIZATION = "none"
try:
    WPS_BIGWIG_NORMALIZED = strtobool(config["wps"]["bigwig_normalized"])
except:
    WPS_BIGWIG_NORMALIZED = False

# defaults : no change in WPS output file suffix
WPS_OUT_SUFFIX = f"_w{WPS_WINDOW}_frg{WPS_MIN}-{WPS_MAX}" if WPS_WINDOW != "120" or WPS_MIN != "120" or WPS_MAX != "180" else ""

SIZE_SELECTION_SUB_FOLDER = "szsel_"+SIZE_SELECTION_MIN+"_"+SIZE_SELECTION_MAX+"/" if SIZE_SELECTION_MIN or SIZE_SELECTION_MAX else ""

# print(f'SIZE_SELECTION_SUB_FOLDER="{SIZE_SELECTION_SUB_FOLDER}"')

try:
    GENOME = config["genome_name"]
except: # default
    GENOME = "hg38"

try:
    CHROM_SIZES = config["chrom_sizes"]
except: # default
    CHROM_SIZES = SCRIPTS+'/data/hg38.chrom.sizes' if GENOME == "hg38" else SCRIPTS+'/data/hg19.chrom.sizes'

try:
    CNA_DO_CONTROL_REFS = config["copy_number_aberration"]["do_control_refs"]
except:
    CNA_DO_CONTROL_REFS = "TRUE"

# still hard coded
RC_NORMALIZATION = "TRUE"
hash_it = False
MIXED_MODEL_FILE_SUFF="mixmod"
HGMM_FILE_SUFF="hgmm"
KEEP_BEDGRAPH = False

CHROMS = config["chromosomes"].split(',')
if not CHROMS: # default: all human chromosome, also for backwar compatibility
    CHROMS = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13",
              "chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM"]

FFT_BINSIZE = str(config["fft"]["binsize"]) if config["fft"]["binsize"] else "10000"
FFT_FREQ_FROM = str(config["fft"]["freq_from"]) if config["fft"]["freq_from"] else "120"
FFT_FREQ_TO = str(config["fft"]["freq_to"]) if config["fft"]["freq_to"] else "280"
try:
    FFT_FREQ_SEL = str(config["fft"]["freq_sel"])
except:
    FFT_FREQ_SEL = "193,196,199"

# import configure as conf
#
# rule configure:
#     output: config['outputDir']+"/"+config['projectName']+'_'+str(config['batch'])+".yaml"
#     params:
#         user_config='/configs/settings.yaml'
#     run:
#         conf.merge_config_and_settings()
