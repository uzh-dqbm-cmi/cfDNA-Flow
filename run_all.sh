#!/bin/bash

for config in configs/all/*.yaml; do
    python cfDNA.py -a do_preprocess -c "$config" -si /cluster/dataset/medinfmk/singularity_images -j16
done
