#!/bin/bash -e 

readonly OUTPUT_DIR="ex2_out"
BINARY=$1

for M in 512 1024
do
    for EPS in "0.01" "0.001"
    do
        for PROC in 32 64 128 256 512
        do
            FILENAME=$1"_"$M"_"$EPS"_"$PROC
            mpisubmit.bg -m vn -w 00:30:00 -n $PROC -stdout $OUTPUT_DIR/$FILENAME.out -stderr $OUTPUT_DIR/$FILENAME.err $BINARY $M $EPS
        done
    done
done
