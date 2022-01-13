#!/bin/bash

zcat Whole_Blood.allpairs.sub.txt.gz | \
    cut -f2-5,7 | \
    sort -k1,1 -k2,2n | \
    uniq | bgzip > Whole_Blood.allpairs.sub.frq.txt.gz
