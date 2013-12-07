#!/bin/bash -e

# Compile and run jobs
# Use: start <source file>
BINARY= echo $1 | sed  -r 's/(.*)\..*/\1/g'
./compile_2ex.sh $1 $BINARY
./run_2ex.sh $BINARY


