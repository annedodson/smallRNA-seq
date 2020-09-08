#!/bin/bash/

# This script creates a summary of the alignment results for all samples.
# USAGE: sh alignment_stats.sh <full pathname of the log directory>

for log in ${1}/*.txt
do
	sample=$(basename $log .txt)
	echo $sample >> ${1}/summary.txt
	cat $log >> ${1}/summary.txt
	echo "" >> ${1}/summary.txt
done
