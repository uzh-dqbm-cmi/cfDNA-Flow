#!/bin/bash

## copy this file as LIQUORICE_singularity.sh, update it according to your environment,
## and add the following line, with the correct 'path' to your .bashrc , or .zshrc :
# 'source /path/LIQUORICE_singularity.sh'
## after that if you want to use it in the same terminal session you need to call:
## 'source ~/.bashrc' or 'source ~/.zshrc' , depending on you shell

# example for Singularity
# use the following command to convert the 'liquorice' Docker to a Singularity image
# 'singularity pull liquorice.sif docker://peneder/liquorice'
#export LIQUORICE_IMG=/path/liquorice.sif
## example for Docker, after 'docker pull peneder/liquorice' :
export LIQUORICE_IMG=peneder/liquorice

function LIQUORICE(){
  # if needed:
#   singularity exec -H $HOME -B $PWD $LIQUORICE_IMG LIQUORICE "$@"
   docker run --rm -v $PWD:$PWD $LIQUORICE_IMG LIQUORICE "$@"
}

function LIQUORICE_summary(){
#  singularity exec -H $HOME -B $PWD $LIQUORICE_IMG LIQUORICE_summary "$@"
  docker run --rm -v $PWD:$PWD $LIQUORICE_IMG LIQUORICE_summary "$@"
}
