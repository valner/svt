#!/bin/bash -e

SOURCE=$1
BINARY=$2

mpixlc_r $SOURCE -O3 -o $BINARY
