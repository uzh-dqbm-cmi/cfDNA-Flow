#!/bin/bash

git clone https://github.com/orisenbazuru/histogram_gmm.git -b fixed_mean && pip install -e ./histogram_gmm

git clone https://github.com/shahcompbio/hmmcopy_utils.git && cd ./hmmcopy_utils && cmake . && make && cd ..

echo "Add this line to your ~/.bashrc|.zshrc:"
echo "export PATH=\$PATH:$PWD/hmmcopy_utils/bin"
echo "To activate it and continue in already opened terminal run (not necessary in the next log-in):"
echo "source ~/.bashrc|.zshrc"
