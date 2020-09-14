#!/bin/bash/

# This script creates a summary of the counts results for all samples.
# USAGE: sh counts_stats.sh <full pathname of the log directory>

for log in ${1}/*.summary
do
	cat $log >> ${1}/summary.txt
	echo "" >> ${1}/summary.txt
done
