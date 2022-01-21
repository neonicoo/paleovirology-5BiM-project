#!/bin/bash -i

# Authors: Adrian Zurmely - Maelle Broustal - Nicolas Mendiboure INSA Lyon 5BIM

##################################### Initial settings  ######################################
##############################################################################################

# donne le repertoire du script : 
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

if [ -d ${SCRIPT_DIR}/plant_239_U100 ] 
then
    :
else
    mkdir plant_239_U100
	cd plant_239_U100

	wget bioinfo.bti.cornell.edu/ftp/program/VirusDetect/virus_database/v239/U100/plant_239_U100.tar.gz
	tar -xzvf plant_239_U100.tar.gz --strip-components=1

	wget bioinfo.bti.cornell.edu/ftp/program/VirusDetect/virus_database/v239/vrl_genbank.info.gz
	tar -xzvf vrl_genbank.info.gz --strip-components=1
	wget bioinfo.bti.cornell.edu/ftp/program/VirusDetect/virus_database/v239/vrl_idmapping.gz
	tar -xzvf vrl_idmapping.gz --strip-components=1
	rm -rf plant_239_U100.tar.gz  vrl_genbank.info.gzvrl_idmapping.gz
	cd ..
fi


if conda env list | grep -q paleogenomic
then
   conda activate paleogenomic
   echo "paleogenomic conda env ON"
else 
	echo "Installing conda env" 
	conda env create --file paleogenomic.yml
	conda clean -a
	conda activate paleogenomic
	echo "paleogenomic conda env ON"
fi


###################################### Pre - processing  #####################################
##############################################################################################

# dossiers : 
# data --> siRNA/ et VANA/
# dans siRNA et VANA --> raw/ (données brutes fastq), QC/, log/ et trimmed_cutadapt/
# QC : multiqc/ multiqc_report.html et fastqc/
# QC/fastqc --> même nom.fastq.html : avant trimming ; _trimmed.html : après trimming 
 

if [ -d ${SCRIPT_DIR}/data/siRNA ] 
then
    :
else
    echo "${SCRIPT_DIR}/data/siRNA doesn't exist"
    exit
fi

if [ -d ${SCRIPT_DIR}/data/VANA ] 
then
    :
else
    echo "${SCRIPT_DIR}/data/VANA doesn't exist"
    exit
fi


printf "\nHello, this is the preprocessing step for VANA and siRNA data.\n\n \nPlease make sure that :\nYour data are in the same root folder as this script meaning :\n"
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

#####################
# siRNA #############
#####################

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

#####################
# VANA ##############
#####################

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

###################################### De novo assembly ######################################
##############################################################################################

function spades_merge ()
{
	# Donne le dossier en argument 1, le fichier en 2e,  \
	# la taille de kmer minimal en 3e argument, \
	# puis la taille de kmer max en 4e argument, \
 	# puis le pas en 5e argument espacés par des espaces.

	if [ $# -eq 0 ]; then
	    echo "Veuillez spécifier le dossier en premier argument,  le fichier en 2e argument, la taille de kmer minimal en 3e argument, puis la taille de kmer max en 4e argument puis le pas en 5e argument espacés par des espaces."
	    exit 1
	fi

	# nom du fichier traité spécifié premier arg:
	FILE=$1
	# dossier spécifié deuxième arg : 
	DOSSIER=$2 
	# valeurs de 11 à 37 par pas de 2
	DEB=$3
	FIN=$4
	PAS=$5

	echo "spades_merge_siRNA.sh en cours"
	echo " dossier $DOSSIER sortie spécifié"
	echo " fichier $FILE entree spécifié"

	echo "kmers de $DEB à $FIN par pas de $PAS"
	cd $DOSSIER

	echo "Assembling $FILE file..."
	for ((i=$DEB;i<$FIN;i+=$PAS)); do 
		echo "spades.py -o contigs_rnaviral/ -s ../../$FILE.fastq -k ${i} --rnaviral "
		spades.py -o contigs_rnaviral/ -s $FILE -k ${i} --rnaviral 
		# récupère le contigs.fasta
		cat contigs_rnaviral/contigs.fasta >> merged_SPAdes.fasta
	done
}

function compress_myfasta ()
{
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
}

#####################
# VANA ##############
##################### 

#spades.py -o outMETA/ -1 Forward_R1.fastq -2 Reverse_R1.fastq --meta

#####################
# siRNA #############
#####################

#spades.py -o test_kmer/ -s ../../190725_SNK268_B_L006_AFVV-15_R1_trimmed.fastq --rnaviral
#ou sans --rnaviral

# test 
# spades.py -o test_kmer/ -s ../../190725_SNK268_B_L006_AFVV-15_R1_trimmed.fastq -k 11,13,15,17,19,21,23,25,27,29,31,33,35,37 --rnaviral
 
# donne le repertoire du script : 
# SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# donne le repertoire du script : 

SCRIPT_DIR=~/test_data/
printf "\nHello, this is the assembly script for VANA and siRNA data.\n\n Please make sure that you ran the preprocessing script : preprocess.sh \n"
UserChoice=0
while [[ $UserChoice != [123] ]]
do
  echo "---------------------------------"
  printf "Which data do you want to assemble ? \n Please type :\n \"1\" for VANA\n \"2\" for siRNA\n \"3\" for both\n \"q\" to quit.\n"
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



#####################
# VANA ##############
#####################

# pour tous les échantillons  VANA/trimmed_cutadapt/
	# bash bbmap/repair.sh in=../BO1-3_Unmapped_R1.fastq_BO1-3_F01.fastq in2=../BO1-3_Unmapped_R2.fastq_BO1-3_F01.fastq
	# spades.py -o output/ -1 Forward_R1.fastq -2 Reverse_R1.fastq
	

	
if [ "$VANA" = true ]; then
	echo "Assembling VANA data"
	FILES_FA="$SCRIPT_DIR/VANA/trimmed_cutadapt/"
	mkdir $FILES_FA/SPAdes	-p
	R1="R1"
	R2="R2"
	for namefileR1 in $FILES_FA/*R1*.fastq
	do
	  #echo "$(basename $namefileR1 .fastq)"
	  namefileR1=$(basename $namefileR1 .fastq)
	  #on remplace le R1 par le R2 :
	  namefileR2=${namefileR1/R1/R2}
	  #echo "R2 = ${namefileR2/R1/R2}"
	  cd $FILES_FA/
	  echo "Creating $namefileR1 directory"
	  mkdir $namefileR1 -p
	  cd $namefileR1
	  echo "Assembling $namefileR1 and $namefileR2 pair-end files..."
	  
	  # on fait le bbmap/repair
	  bash ~/scripts/bbmap/repair.sh in=$FILES_FA/$namefileR1.fastq in2=$FILES_FA/$namefileR2.fastq out1=$FILES_FA/$namefileR1/$namefileR1.repaired_R1.fastq out2=$FILES_FA/$namefileR1/$namefileR2.repaired_R2.fastq outs=$FILES_FA/$namefileR1/$namefileR1.unpaired.fastq
	  # lance spades : 
	  spades.py -o outputMETA/ -1 $FILES_FA/$namefileR1/$namefileR1.repaired_R1.fastq -2 $FILES_FA/$namefileR1/$namefileR2.repaired_R2.fastq -s $FILES_FA/$namefileR1/$namefileR1.unpaired.fastq --meta
	  
	  #pour idba_ud on merge 1 et 2 en un seul fichier sans les unpaired :
	  fq2fa --merge $FILES_FA/$namefileR1/$namefileR1.repaired_R1.fastq $FILES_FA/$namefileR1/$namefileR2.repaired_R2.fastq $FILES_FA/$namefileR1/$namefileR1.mergedR1-R2.fastq
	  idba_ud -r $FILES_FA/$namefileR1/$namefileR1.mergedR1-R2.fastq --num_threads 8 -o idba_ud_out

	  done
fi
	

#####################
# siRNA #############
#####################

#pour tous les échantillons  siRNA/trimmed_cutadapt/
	# pour tous les fichiers test_data/siRNA/trimmed_cutadapt/SPAdes/BO1_15/test_kmer/K*/final_contigs.fasta :
	# echo le fichier >> merged_fasta_kmer.fasta
	
	
if [ "$siRNA" = true ]; then
	echo "Assembling siRNA data"
	FILES_FA="$SCRIPT_DIR/siRNA/trimmed_cutadapt/"
	for f in $FILES_FA*.fastq
	do
	  cd $FILES_FA/
	  f=$(basename $f .fastq)
	  echo "Creating $f directory"
	  mkdir $f -p
	  cd $f
	  mkdir SPAdes -p
	  # lance spades_merge_siRNA sur fichier fasta 1e arg, dossier sortie 2e arg, kmers de 11 à 37 pas de 2 :
	  spades_merge $FILES_FA/$f.fastq $FILES_FA/$f/SPAdes/ 11 37 2
	  #bash idba_merge_siRNA.sh . $f 11 37 2
	  
	  # compress les fasta avec le script de denis si y a redondance : 
	  cd $FILES_FA/$f/SPAdes/
      cp merged_SPAdes.fasta compressed_SPAdes.fasta
      compress_myfasta compressed_SPAdes.fasta
	done
fi

########################################## Blast #############################################
##############################################################################################






FILEPATH=$2 #contigs to blast
FILE="${FILEPATH##*/}"
FILENAME="${FILE%.*}"

function myblastn ()
{
	blastn -query $1 \
		  	-db $2 \
		  	-out $3\
		  	-num_threads 8 \
		  	-evalue 0.01 \
		  	-outfmt 6
}

function myblastx ()
{
	blastx -query $1 \
			-db $2\
			-out $3\
			-num_threads 8 \
			-evalue 0.01 \
			-outfmt 6
}


while [ -n "$1" ]; do # while loop starts

	case "$1" in

	-virusdetect) 
		echo "You're going to blast{n,x} your contigs over the virus plants database"
		echo -n "Continue? [y/n] " 
		read doit
			case $doit in 
				y|Y) 
					# $3 : virusdetect_db

					myblastn $FILEPATH $3vrl_Plants_239_U100 ${FILENAME}_virusdetect_blastn.txt
					echo "#### blastN on virusdetect db - success ####"

					myblastx $FILEPATH $3vrl_Plants_239_U100_prot ${FILENAME}_virusdetect_blastx.txt
					echo "#### blastX on virusdetect db- success ####"

					echo "#### Virus identification ####"

					if [[ -s ${FILENAME}_virusdetect_blastn.txt ]]
					then 
						python3 ./blastn_virus_identity.py ${FILENAME}_virusdetect_blastn.txt $3 ${FILENAME}_virusdetect_blastn_taxon
					fi

					if [[ -s ${FILENAME}_virusdetect_blastx.txt ]]
					then 
						python3 ./blastx_virus_identify.py ${FILENAME}_virusdetect_blastx.txt $3 ${FILENAME}_virusdetect_blastx_taxon
					fi

					echo "#### DONE ####" ;;

				n|N) echo "Cancel ";; 
				*) echo "Cancel" ;; 
			esac
			shift
		;;

	-nrnt)
		echo "You're going to blast{n,x} your contigs over the nr/nt database"
		echo -n "Continue? [y/n] " 
		read doit2
			case $doit2 in 
				y|Y) 
					# $3 : nr db
					# $4 : nt db

					myblastn $FILEPATH $3 ${FILENAME}_nr_blastn.txt
					echo "#### blastN on nr db - success ####"

					myblastx $FILEPATH $3 ${FILENAME}_nt_blastx.txt
					echo "#### blastX nr db - success ####"

					myblastn $FILEPATH $4 ${FILENAME}_nt_blastn.txt
					echo "#### blastN on nt db - success ####"

					myblastx $FILEPATH $4 ${FILENAME}_nt_blastx.txt
					echo "#### blastX nt db - success ####"

					echo "#### DONE ####" ;;

				n|N) echo "Cancel ";; 
				*) echo "Cancel" ;;
			esac
			shift
		;;
	esac
	shift

done

conda deactivate
echo "paleogenomic conda env OFF"