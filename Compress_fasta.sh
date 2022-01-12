#!/bin/bash

# Authors: Linux de Folis / Denis Filloux - CIRAD - UMR PHIM - 2021

# Usage: ./Compress_fasta.sh File.fa
# Remove duplicated sequences in a multifasta file
# Be careful: sequence IDs must be different
# Output: $1_Compressed.fa

module load fastx_toolkit/0.0.14
module load vsearch/2.14.0

echo "./Compress_fasta.sh "$1
echo $(grep -c "^>" $1)" sequences in "$1

vsearch --cluster_fast $1 --centroids $1_Compressed.fa --iddef 0 --id 1.00 --strand both --qmask none --fasta_width 0 --minseqlength 1 --maxaccept 0 --maxreject 0

fasta_formatter -w 0 -i $1 -o $1.tab -t
fasta_formatter -w 0 -i $1_Compressed.fa -o $1_Compressed.fa.tab -t
cut -f2 $1_Compressed.fa.tab | rev | tr atgcATGC tacgTACG > $1_Compressed.fa.tab.tab
paste $1_Compressed.fa.tab $1_Compressed.fa.tab.tab > $1_Compressed.fa.tab.tab.tab
rm $1_Compressed.fa.tab
rm $1_Compressed.fa.tab.tab
awk -F $"\t" '{print ">"$1"|RC\n"$3}' $1_Compressed.fa.tab.tab.tab > $1_Compressed.fa.rc 
rm $1_Compressed.fa.tab.tab.tab
cat $1_Compressed.fa $1_Compressed.fa.rc > $1_Compressed.fa.all
rm $1_Compressed.fa.rc
old_IFS=$IFS
IFS=$'\t'
> $1_Compressed.fa.tab
while read c1 c2
	do
	printf "$c1\t$c2\t" >> $1_Compressed.fa.tab
	grep -c $c2 $1_Compressed.fa.all >> $1_Compressed.fa.tab
done < $1.tab
IFS=$old_IFS
rm $1.tab
rm $1_Compressed.fa.all
awk -F $"\t" '{if ($3==0) print ">"$1"\n"$2}' $1_Compressed.fa.tab > $1_Compressed.fa.more
rm $1_Compressed.fa.tab
echo $(grep -c "^>" $1_Compressed.fa.more)" readded sequences"
cat $1_Compressed.fa.more >> $1_Compressed.fa
rm $1_Compressed.fa.more 

fasta_formatter -w 0 -i $1_Compressed.fa -o $1_Compressed.fa.tab -t 
awk -F $"\t" '{print $1"\t"$2"\t"length($2)}' $1_Compressed.fa.tab | sort -t $'\t' -n -r -k3,3 | awk -F $"\t" '{print ">"$1"\n"$2}' > $1_Compressed.fa
rm $1_Compressed.fa.tab

echo $(grep -c "^>" $1_Compressed.fa)" sequences in "$1_Compressed.fa

