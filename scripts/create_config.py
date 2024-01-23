import argparse
import pandas as pd
from pathlib import Path
import yaml


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-fq_dir", dest="fastq_dir", required=True)
    parser.add_argument("-meta", dest="meta_csv",  required=True)
    parser.add_argument("-ref", dest="genome", required=True)
    parser.add_argument("-s", dest="settings_file", default="settings.yaml", required=False)
    parser.add_argument("-fwd", dest="forward", default="_R1", required=False)
    parser.add_argument("-rvr", dest="reverse", default="_R2", required=False)
    parser.add_argument("-n", dest="projectName",  default="Test", required=False)
    parser.add_argument("-o", dest="outputDir", default=".", required=False)
    parser.add_argument("-cond", dest="condition", default="tissue", required=False)
    parser.add_argument("-C", dest="capture", default="", required=False)
    parser.add_argument("-B", dest="bait", default="", required=False)
    parser.add_argument("-dbSNP", dest="dbSNP", default="", required=False)
    parser.add_argument("-gnomad", dest="gnomad", default="", required=False)
    parser.add_argument("-pon", dest="PON", action="store_true", required=False)
    parser.add_argument("-germRes", dest="germRes", action="store_true", required=False)
    parser.add_argument("-L", dest="intList", default="", required=False)
    return parser


def process_meta_file(meta_csv, condition='tissue'):
    df = pd.read_csv(meta_csv).applymap(str)
    if condition not in df.columns:
        df[condition] = 'normal'
    return df[["SampleID", "PatientID", condition]]


def create_sample_dict(df, fastq_dir, fwd, rvr):
    files = [str(f) for f in Path(fastq_dir).iterdir()]
    print(df)
    sample_dict = df.set_index("SampleID").T.to_dict()
    for sample in sample_dict.keys():
        for file in files:
            if str(sample) in file:
                if fwd in file:
                    sample_dict[sample]["forward"] = file
                elif rvr in file:
                    sample_dict[sample]["reverse"] = file
                else:
                    print("Cannot find fastq file for sample")
    return sample_dict


def create_patient_dict(sample_dict, condition='tissue'):
    patient_dict = {}
    for k, v in sample_dict.items():
        pid = v["PatientID"]
        if pid not in patient_dict.keys():
            patient_dict[pid] = {}
        patient_dict[pid][v[condition]] = k
    return patient_dict


def read_settings(settings_file="settings.yaml"):
    with open(settings_file, 'r') as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            return exc

def create_config_file(parser):
    args = parser.parse_args()  # todo create ouput directory if doesn't exist
    df = process_meta_file(args.meta_csv, args.condition)
    sample_dict = create_sample_dict(df, args.fastq_dir, args.forward, args.reverse)
    patient_dict = create_patient_dict(sample_dict, args.condition)
    settings = read_settings(args.settings_file)
    final_dict = {**{"projectName": args.projectName,
                     "outputDir": args.outputDir,
                     "genome": args.genome,
                     "capture": args.capture,
                     "bait": args.bait,
                     "dbSNP": args.dbSNP,
                     "gnomad": args.gnomad,
                     "PON": args.PON,
                     "withGermRes": args.germRes,
                     "samples": sample_dict,
                     "patients": patient_dict,
                     "intList": args.intList,
                     }, **settings}
    with open(str(Path(args.outputDir)/"griffin_config.yaml"), 'w') as outfile:
        yaml.dump(final_dict, outfile, default_flow_style=False)
    return final_dict


if __name__ == "__main__":
    create_config_file(get_args())
