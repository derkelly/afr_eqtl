#!/bin/bash

PED=$1
K=$2

~/bin/admixture --cv -s 4690 $PED $K
