#!/bin/bash

A=$1
B=$2
OUT=$3

bedtools intersect -a $A -b $B | bgzip > $OUT
