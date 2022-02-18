#!/bin/bash

if [[ $# -ne 3 ]]; then
    echo "ERROR: script accepts 3 arguments: particles, destination, cores">&2
    exit 2
fi

cd build

DEST=../output/$2_scaling.out

export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export OMP_NUM_THREADS=$3

echo $3 cores >> $DEST
./openmp -s 1 -n $1 >> $DEST
