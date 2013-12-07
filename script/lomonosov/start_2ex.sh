#!/bin/bash -e

# Compile and run jobs
# Use: start <source file>
module add openmpi/1.5.5-icc

module add intel/13.1.0

module add mkl/4.0.2.146

./compile_2ex_enter.sh 
./run_2ex.sh "sor"


