#!/bin/bash

cd /data/auyar/ATAC-Seq_Data/Stitzel/Alzheimer_ATACseq_human

# files=`find ./ -name *_sorted_peaks.bed`
files=`find ./ -name *_sorted_peaks.broadPeak`

for idx in $files; do
# cp $idx /data/xwang/ATAC
  echo $idx
done
