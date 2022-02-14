#!/bin/bash

DATE=$(date +'%m-%d_%H-%M')
DIRNAME="output/$DATE"
if [[ $# -eq 1 ]]; then
    DIRNAME="${DIRNAME}_${1}"
    echo $DIRNAME
    mkdir $DIRNAME
fi

# cd build
for f in build/job*; do
    # echo $f
    sbatch $f
done
