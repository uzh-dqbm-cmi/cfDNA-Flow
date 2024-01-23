#!/bin/bash


if samtools --version | grep -c "no version information available"
then
  echo "ERROR: samtools install failed"
  exit 1
else
  echo "samtools install OK"
fi


samtools --version | if grep -q "Samtools compilation details" ; then if grep -q "no version information available" ; then echo "ERROR: samtools install failed" && exit 1; else echo "samtools install OK"; fi; else echo "ERROR: samtools install failed" && exit 1; fi
skewer --help | if  grep -q "USAGE:" ; then echo "skewer install OK"; else echo "ERROR: skewer install failed" && exit 1; fi
fastqc --help | if  grep -q "A high throughput sequence QC analysis tool" ; then echo "fastqc install OK"; else echo "ERROR: fastqc install failed" && exit 1; fi
multiqc --help | if  grep -q "Usage:" ; then echo "multiqc install OK"; else echo "ERROR: multiqc install failed" && exit 1; fi
bedtools --help | if  grep -q "Usage:" ; then echo "bedtools install OK"; else echo "ERROR: bedtools install failed" && exit 1; fi
deeptools --help | if  grep -q "usage:" ; then echo "deeptools install OK"; else echo "ERROR: deeptools install failed" && exit 1; fi
# bedGraphToBigWig | if  grep -q "usage:" ; then echo "bedGraphToBigWig install OK"; else echo "ERROR: bedGraphToBigWig install failed" && exit 1; fi
gatk --help | if  grep -q "Usage" ; then echo "gatk install OK"; else echo "ERROR: gatk install failed" && exit 1; fi
bismark --help | if  grep -q "USAGE:" ; then echo "bismark install OK"; else echo "ERROR: bismark install failed" && exit 1; fi
bedops --help | if  grep -q "USAGE:" ; then echo "bedops install OK"; else echo "ERROR: bedops install failed" && exit 1; fi
bowtie2 --help | if  grep -q "Usage:" ; then echo "bowtie2 install OK"; else echo "ERROR: bowtie2 install failed" && exit 1; fi
bwa | if egrep "Usage:" ; then echo "bwa install OK"; else echo "ERROR: bwa install failed" && exit 1; fi
# singularity shell -H $HOME ./cfDNA_pipeline.sif
# pytest -q --inst_folder / /test_cfDNA.py
# snakemake -s Snakefile --configfile test/test_cfDNA_pipeline.yaml -j 2 do_preprocess -np
# /test_cfDNA.py  # cannot create results folder
# pytest -q --inst_folder / /test_cfDNA.py::test_smk_do_preprocess
