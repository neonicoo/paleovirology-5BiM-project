#!/bin/bash

#preprocessing.sh

# Authors: Adrian Zurmely - Maelle Broustal - Nicolas Mendiboure INSA Lyon 5BIM

# Usage: ./preprocessing.sh
# Output: trimmed data in VANA/trimmed_cutadapt and siRNA/trimmed_cutadapt MultiQC in VANA/QC and siRNA/QC

# dossiers : 
# data --> siRNA/ et VANA/
# dans siRNA et VANA --> raw/ (données brutes fastq), QC/, log/ et trimmed_cutadapt/
# QC : multiqc/ multiqc_report.html et fastqc/
# QC/fastqc --> même nom.fastq.html : avant trimming ; _trimmed.html : après trimming 
 

# donne le repertoire du script : 
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )



printf "\nHello, this is the preprocessing script for VANA and siRNA data.\n\n \nPlease make sure that :\nYour data are in the same root folder as this script meaning :\n"
printf "in the data/ folder there should be :\n"
printf " -this script \"preprocessing.sh\" \n"
printf " -VANA/raw folders with raw VANA data .fastq R1 and R2 paired end files \n"
printf " -siRNA/raw folders with raw siRNA data .fastq files \n"
UserChoice=0
while [[ $UserChoice != [123] ]]
do
  echo "---------------------------------"
  printf "Which data do you want to preprocess ? \n Please type :\n \"1\" for VANA\n \"2\" for siRNA\n \"3\" for both\n \"q\" to quit.\n"
  read -p 'Data to preprocess : ' UserChoice
  if [[ $UserChoice == [qQ] ]]; then
    break
  fi
done


siRNA=false
VANA=false

case  $UserChoice in 

	3)
	  siRNA=true
	  VANA=true
	;;
	2)
	  siRNA=true
	;;
	1)
	  VANA=true
	;;
esac

#######
#siRNA

if [ "$siRNA" = true ]; then
	echo "Preprocessing siRNA data"
	cd $SCRIPT_DIR
	cd siRNA/
	FILES="raw/"
	mkdir -p QC
	mkdir -p QC/fastqc
	mkdir -p log
	mkdir -p trimmed_cutadapt
	#fastqc 
	fastqc -t 6 raw/*.fastq -o QC/fastqc > log/fastqc_pre.txt

	#cutadapt
	ADAPT_ILLUMINA_siRNA=TGGAATTCTC
	ADAPT_ILLUMINA_siRNA_RC=CACCCGAGAATTCCA

	for f in $FILES*.fastq
	do
	  echo "$(basename $f .fastq)"
	  f=$(basename $f .fastq)
	  echo "Processing $f file..."
	  # take action on each file. $f store current file name
	  cutadapt -m 10 -q 30 -a $ADAPT_ILLUMINA_siRNA -a $ADAPT_ILLUMINA_siRNA_RC -o trimmed_cutadapt/${f}_trimmed.fastq  $FILES/$f.fastq > log/cutadapt_$f.txt
	done

	#fastqc
	fastqc -t 6 trimmed_cutadapt/*_trimmed.fastq -o QC/fastqc > log/fastqc_post.txt

	#multiqc :
	multiqc -s -f QC/fastqc -o QC > log/multiqc.txt
fi

#######
#VANA :
if [ "$VANA" = true ]; then
	cd $SCRIPT_DIR
	echo "Preprocessing VANA data"
	cd VANA/
	FILES="raw/"
	mkdir -p QC
	mkdir -p QC/fastqc
	mkdir -p log
	mkdir -p trimmed_cutadapt
	#rearranging filenames : ".fastq.gz" in a middle of a filename fastqc report can mess up with multiqc analysis
        cd  $FILES
        for file in *
        do
	  if [[ "$file" == *".fastq.gz"* ]];then
            mv "$file" "${file/.fastq.gz/}"
	  fi
	done
        cd $SCRIPT_DIR/VANA/

	#fastqc
	fastqc -t 6 raw/*.fastq -o QC/fastqc > log/fastqc_pre.txt

	#cutadapt
	ADAPT_R1="AGATCGGAAGAGCAC"
	ADAPT_R2="AGATCGGAAGAGCGT"
	ADAPT_RC_R1="GTGCTCTTCCGATCT"
	ADAPT_RC_R2="ACGCTCTTCCGATCT"

	R1="R1"
	R2="R2"
	for namefileR1 in $FILES*R1*.fastq
	do
	  #echo "$(basename $namefileR1 .fastq)"
	  namefileR1=$(basename $namefileR1 .fastq)
	  #on remplace le R1 par le R2 :
	  namefileR2=${namefileR1/R1/R2}
	  #echo "R2 = ${namefileR2/R1/R2}"
	  echo "Processing $namefileR1 and $namefileR2 pair-end files..."
	  cutadapt -m 20 -q 30 -u 24 -U 24 -a $ADAPT_R1 -a $ADAPT_RC_R1 -A $ADAPT_R2 -A $ADAPT_RC_R2 -o trimmed_cutadapt/${namefileR1}_trimmed.fastq -p trimmed_cutadapt/${namefileR2}_trimmed.fastq $FILES/$namefileR1.fastq $FILES/$namefileR2.fastq > log/cutadapt_$namefileR1.txt
	done
	#fastqc
	fastqc -t 6 trimmed_cutadapt/*_trimmed.fastq -o QC/fastqc > log/fastqc_post.txt
	#multiqc :
	multiqc -s -f QC/fastqc/*zip -o QC > log/multiqc.txt

fi
