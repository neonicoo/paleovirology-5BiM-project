#!/bin/bash

# Authors: Linux de Folis / Denis Filloux - CIRAD - UMR PHIM - 2018

# Usage : ./Redondant_fastq.sh File.fastq
# Remove all redondant sequences (based on the sequence and quality) in fastq file
# Output: Non_redondant_$1

echo "./Redondant_fastq.sh" $1

echo $(expr $(wc -l $1 | cut -d " " -f1) / 4 )" sequences in "$1
cat $1 | paste - - - - | awk -F $"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$2$4}' > $1.tab
sort -t $'\t' -k5,5 -u $1.tab > $1.sort.tab
rm $1.tab
awk -F $"\t" '{print $1"\n"$2"\n"$3"\n"$4}' $1.sort.tab > Non_redondant_$1
rm $1.sort.tab

echo $(expr $(wc -l Non_redondant_$1 | cut -d " " -f1) / 4 )" sequences in "Non_redondant_$1
