#!/bin/bash

readonly OUTPUT_DIR="2ex_out"
BINARY=$1
if [ ! -d $OUTPUD_DIR ];
then
mkdir $OUTPUT_DIR
fi
for M in 512 1024
do
    for EPS in "0.01" "0.001"
    do
        for PROC in 32 64 128 256
        do
            FILENAME=$1"_"$M"_"$EPS"_"$PROC
            sbatch -n$PROC -o $OUTPUT_DIR/$FILENAME.out -e $OUTPUT_DIR/$FILENAME.err ompi $BINARY $M $EPS
        done
    done
done
