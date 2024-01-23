#! /usr/bin/env python

import pandas as pd
import os
import re
import sys
import glob
# Helper script to extract sample infos from meta data for the public data.
#
# Study name: Genome wide cell-free DNA fragmentation in patients with cancer
#    Study accession: EGAS00001003611
#    Dataset accession: EGAD00001005339

sample_sile_file = 'Sample_File.map.tsv'
meta_info_file = 'Run_Sample_meta_info.map.csv'
samples_folder = '.'

output_info_file = 'samples_info.tsv'


sample_sile_df = pd.read_table(sample_sile_file, names=['id', 'sample_dir'], usecols=[0, 3])
meta_info_df = pd.read_csv(meta_info_file, delimiter=';', names=['subject_id', 'gender', 'phenotype'], usecols=[0, 1, 2])

with open(output_info_file, 'w') as f:
    # print("sample_dir\tsubject_id\tgender\tphenotype")
    f.write("sample_dir\tsubject_id\tgender\tphenotype\n")
    for rootdir, dirs, files in os.walk(samples_folder):
        for sample_dir in dirs:
            bam = glob.glob(sample_dir + "/*.bam")
            if not bam:
                continue
            id_df = sample_sile_df[sample_sile_df['sample_dir'] == sample_dir]
            if not len(id_df):
                sys.stderr.write(f"{sample_dir}\tmissing information\n")
                continue
            id = sample_sile_df[sample_sile_df['sample_dir'] == sample_dir].iloc[0]['id']
            cleaned_id = re.sub('_(.+)', '', id)
            subject_id_pref = f'subject_id={cleaned_id}_'
            # print(f"sample_dir:{sample_dir}\tid: {id}\tsubject_id_pref: {subject_id_pref}")
            sample_meta = meta_info_df[meta_info_df['subject_id'].str.contains(subject_id_pref)].iloc[0]
            # print(f"{sample_dir}\t{sample_meta['subject_id'].split('=')[1]}\t{sample_meta['gender'].split('=')[1]}\t{sample_meta['phenotype'].split('=')[1]}")
            f.write(f"{sample_dir}\t{sample_meta['subject_id'].split('=')[1]}\t{sample_meta['gender'].split('=')[1]}\t{sample_meta['phenotype'].split('=')[1]}\n")

