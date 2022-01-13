#!/bin/bash

JUNCS=$1
OUT=$2

python ~/bin/leafcutter/clustering/leafcutter_cluster.py \
-j $JUNCS \
-m 50 \
-o $OUT \
-l 500000
