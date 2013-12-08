#!/bin/bash -e
readonly MAPFILE="map.txt"
BINARY=`echo $1 | sed  -r 's/(.*)\..*/\1/g'`
./compile_2ex.sh $1 $BINARY
./generate_map $MAPFILE
readonly OUTPUT_DIR="ex2_out"
BINARY=$1

for M in 512 1024
do
    for EPS in "0.01" "0.001"
    do
        PROC=512
        FILENAME=$1"_"$M"_"$EPS"_"$PROC_randmap
        mpisubmit.bg -e \"MPIRUN_MAPFILE=$MAPFILE\" -m vn -w 00:30:00 -n $PROC -stdout $OUTPUT_DIR/$FILENAME.out -stderr $OUTPUT_DIR/$FILENAME.err $BINARY $M $EPS
    done
done


