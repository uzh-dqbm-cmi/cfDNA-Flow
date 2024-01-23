include: "config.smk"

import numpy as np
import pandas as pd
import re
from snakemake.remote import FTP
from snakemake.utils import validate

global CONTROL_SAMPLES
global WORKDIR
global hash_it

global samples
global control_samples
global control_RefName
global treated_samples
global treated_samples_names
global treated_samples_RefName

ftp = FTP.RemoteProvider()


samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

validate(samples, schema="schemas/samples.schema.yaml")


import pandas as pd
from hashlib import md5

def aligner(wildcards):
    if config['bwa_algorithm'] == 'mem':
        return expand(WORKDIR + "/BAM/{s}.mem_sorted.bam", s=wildcards.sample)
    elif config['bwa_algorithm'] == 'mem2':
        return expand(WORKDIR + "/BAM/{s}.mem2_sorted.bam", s=wildcards.sample)
    elif config['bwa_algorithm'] == 'aln':
        return expand(WORKDIR + "/BAM/{s}.aln_sorted.bam", s=wildcards.sample)


def gen_ref_name(file_names, prefix):
    ref_str = "_".join(sorted(file_names))
    # print(f'len(file_names)={len(file_names)}')
    ref_id = md5(ref_str.encode()).hexdigest() if hash_it or len(file_names) > 7 else ref_str
    ref_name = make_name(prefix + ref_id)
    print(f"Using refName {ref_name} for {ref_str}")
    return ref_name


def make_name(name: str):
    new_name = re.sub('[^a-zA-Z0-9àèçÀÈÇäöüÄÖÜ_]', '.', name)
    new_name = f'X{name}' if new_name[0].isdigit() or new_name[0] == '_' or new_name[0] == '.'else new_name
    if new_name != name:
        print(f'make_names: name {name} changed to {new_name}')
    return new_name

def make_names(names):
    """
    Implements the same transformation like in R's make.names()
    :param names: list of strings to be made as column conform names
    :return: column conform names: only letters, digits, '.' and '_' allowed. Starting only with a letter, otherwise an 'X' is added as a prefix.
    """
    if isinstance(names, (list, tuple, np.ndarray)):
        new_names = [make_name(name) for name in names]
    else:
        new_names = make_name(names)
    return new_names


control_samples = pd.read_csv(CONTROL_SAMPLES, header=None) #
control_RefName = gen_ref_name(control_samples[0].tolist(), "Ref_")   #"Ref_"+ "_".join(control_samples[0].tolist())

# treated_samples = list(set(samples_df['sample_name'].values) - set(control_samples[0].tolist()))
treated_samples = [x for x in samples['sample_name'].values if x not in control_samples[0].tolist()]
treated_samples_RefName = gen_ref_name(treated_samples, "ALL_")
treated_samples_names = [make_name(x) for x in treated_samples]
