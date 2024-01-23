#!/usr/bin/env python3
# deprecated: instead use configure.py
import argparse
import subprocess
import shlex
import shutil
# from pathlib import Path
import yaml
# from ruamel import yaml
import os
import sys
import pprint

import configure

# Options


def get_args():
    parser = argparse.ArgumentParser("cfDNA Preprocessing Pipeline\n")
    parser.add_argument("-a", "--analysis",
                        help="Options: preprocess, call, annotate, wps, fft", required=True)
    parser.add_argument("-c", "--config", nargs='?', type=str,
                        help="Config file", required=False)
    parser.add_argument('-j', '--j', help="Parallel param in snakemake",
                        type=int, default=1, required=False)
    parser.add_argument('-np', '--np', help="Dry run", action='store_true',
                        required=False)
    parser.add_argument("-l", "--local", help='Run locally or on cluster',
                        action='store_true', required=False)
    parser.add_argument('-u', '--unlock', action='store_true', required=False)
    parser.add_argument('-si', '--singularity_images', default='singularity_images',
                        help='Location of singularity images folder.')
    return parser

def build_image(singularity_images, analysis, local=True):
    analysis_to_image = {'do_preprocess': ['cfDNA_pipeline'],
                         'do_qc': ['cfDNA_pipeline'],
                         'do_bedprocess': ['cfDNA_pipeline'],
                         'do_build': ['cfDNA_pipeline'],
                         'do_wps': ['cfDNA_pipeline']}
    image = analysis_to_image[analysis][0]
    image_file = f"{image}.sif"
    image_path = os.path.join(singularity_images, image_file)
    if local:
        if not os.path.isfile(image_path):
            print(f"No image found, building {image_path} image")
            cmd = shlex.split(
                f"sudo singularity build {image}.def {image_path}")
            subprocess.call(cmd)
        else:
            print(f"{image_path} image found")
        return image_path
    else:
        if analysis == 'build':
            print("Can't build on LeoMed")
            sys.exit()
        else:
            if not os.path.isfile(image_path):
                print(f"No image {image_path} found, Exiting")
                sys.exit()
            else:
                return image_path



def clean(config):
    user_configs = yaml.load(open(config), Loader=yaml.FullLoader)
    print("Removing {} ".format(user_configs['outputDir']))
    shutil.rmtree(user_configs['outputDir'])


def main(args, test=False):

    print(args)
    analysis = args.analysis

    image = build_image(args.singularity_images, analysis, args.local)

    if analysis == 'build':
        return "Done"

    user_config = args.config
    j = args.j
    np = '-np' if args.np else ''
    unlock = '--unlock' if args.unlock else ''
    local = args.local
    new_config = configure.merge_config_and_settings(user_config)
    pp = pprint.PrettyPrinter(indent=2)
    pp.pprint(new_config)

    if local:
        cmd_str = "singularity exec {} snakemake --configfile {} " \
                  "-j {} " \
                  "{} {} {}".format(image, new_config, j, np, analysis, unlock)
        print(cmd_str)
        subprocess.call(shlex.split(cmd_str))
    else:
        new_config_yaml = yaml.load(open(new_config), Loader=yaml.FullLoader)
        leomed = new_config_yaml['leomed']
        nodes = leomed['nodes']
        mem = leomed['mem']
        wallTime = leomed['walltime']
        dataDir = leomed['data']
        # added, because LeoMed does not resolve $HOME
        home = os.environ['HOME']
        jobname = new_config_yaml['projectName']  # 'cfDNA_pipeline'
        cmd_str1 = f"bsub -n {nodes} -R 'rusage[mem={mem}]' -W {wallTime} -J {jobname} -o '{jobname}_%J.log' "
        cmd_str2 = f"singularity exec -H {home} -B {dataDir} {image} "
        # cmd_str2 = "singularity exec -H {} -B {}:/opt/data {} ".format(home, dataDir, image)
        # line below removed because LeoMed does not resolve $HOME
        # cmd_str2 = "singularity exec -H $HOME -B {}:/opt/data {} ".format(dataDir, image)
        # todo make sure local path to image works on leomed
        pwd = os.path.dirname(os.path.realpath(__file__))  # os.getcwd()
        snakemake_file = os.path.join(pwd, "Snakefile")

        cmd_str3 = f"snakemake -s {snakemake_file} --configfile {new_config} -j {j} {np} {analysis} {unlock}"
        cmd_str = cmd_str1 + cmd_str2 + cmd_str3
        print(cmd_str)
        # todo test this on leomed
        if test:
            print(f'run command: {shlex.split(cmd_str)}')
        else:
            subprocess.call(shlex.split(cmd_str))


if __name__ == "__main__":
    args = get_args().parse_args()
    if args.analysis == 'clean':
        clean(args.config)
    else:
        main(args)
