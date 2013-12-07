#!/bin/bash -e

SOURCE=$1
BINARY=$2

mpicc $SOURCE -O3 -o $BINARY
