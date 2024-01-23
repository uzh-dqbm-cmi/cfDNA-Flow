#!/usr/bin/env python3

import yaml
import os
import subprocess
import argparse


def get_args():
    parser = argparse.ArgumentParser("cfDNA Preprocessing Pipeline.")
    parser.add_argument("-c", "--user_config", required=True,
                        help="User config file to be merged with the default settings.yaml config")
    parser.add_argument("-o", "--out_config",  # default='./config.yaml',
                        help="Output merged config file.")
    return parser


def merge_dicts(x, y):
    z = x.copy()
    z.update(y)
    clean_z = none_to_empty_str(z)
    return clean_z


def none_to_empty_str(dict):
    return {k: ('' if v is None else v) for k, v in dict.items()}


def merge_config_and_settings(config, out_config=None):
    user_configs = yaml.load(open(config), Loader=yaml.FullLoader)
    pwd = os.path.dirname(os.path.realpath(__file__))  # os.getcwd()
    settings_yaml = os.path.join(pwd, "configs", "settings.yaml")
    settings_config = yaml.load(open(settings_yaml), Loader=yaml.FullLoader)
    subprocess.call(["mkdir", "-p", user_configs['outputDir']])
    new_config_path = out_config if out_config else os.path.join(user_configs['outputDir'], user_configs['projectName']
                                                                 + '_' + str(user_configs['batch']) + ".yaml")
    if os.path.isfile(new_config_path):
        return new_config_path
    else:
        merged_config = merge_dicts(settings_config, user_configs)  # new way: real merge
        with open(new_config_path, 'w') as outfile:
            yaml.dump(merged_config, outfile, default_flow_style=False)
        #with open(new_config_path, "w") as new_config:
        #    new_config.write(open(settings_yaml, "r").read())
        #    new_config.write(open(config, "r").read())
        return new_config_path


if __name__ == "__main__":
    args = get_args().parse_args()
    merge_config_and_settings(args.user_config)
