#!/bin/bash -e

# warning: filenames mandatory to sor.c
SOURCE=sor.c
BINARY=sor
module add openmpi/1.5.5-icc

module add intel/13.1.0

module add mkl/4.0.2.146
mpicc $SOURCE -O3 -o $BINARY
logout
