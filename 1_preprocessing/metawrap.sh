#!/bin/bash

conda activate metaWRAP

mkdir READ_QC
for F in RAW_READS/*_1.fastq; do 
	R=${F%_*}_2.fastq
	BASE=${F##*/}
	SAMPLE=${BASE%_*}
	metawrap read_qc -1 $F -2 $R -t 6 -o READ_QC/$SAMPLE
done > metawrapOUT.txt 2>&1 # obtains log file
