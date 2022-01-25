#!/bin/bash

# Authors: Adrian Zurmely - Maelle Broustal - Nicolas Mendiboure INSA Lyon 5BIM

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
DATA_DIR="${SCRIPT_DIR}/data"
SIRNA_DIR="$DATA_DIR/siRNA"
VANA_DIR="$DATA_DIR/VANA"
DATABASES_DIR="${SCRIPT_DIR}/databases"


# dossiers : 
# data --> siRNA/ et VANA/
# dans siRNA et VANA --> raw/ (données brutes fastq), QC/, log/ et trimmed_cutadapt/
# QC : multiqc/ multiqc_report.html et fastqc/
# QC/fastqc --> même nom.fastq.html : avant trimming ; _trimmed.html : après trimming 

if [ -d "${SIRNA_DIR}/trimmed_cutadapt" ] 
then
    :
else
    echo "${SIRNA_DIR}/trimmed_cutadapt doesn't exist"
    exit
fi

if [ -d "${VANA_DIR}/trimmed_cutadapt" ] 
then
    :
else
    echo "${VANA_DIR}/trimmed_cutadapt doesn't exist"
    exit
fi

printf "\nThis part allows you to execute a first sorting of your datas: you can choose to eliminate some of your reads to ease and speed up the next steps.\n\n"
printf "\n This step isn't mandatory. Please choose if you want to execute it."
UserChoice=0
while [[ $UserChoice != [yn] ]]
do
  printf "\n ---------------------------------\n"
  printf "Execute kraken reads classification ?"
  read -p '[y/n]' UserChoice
  if [[ $UserChoice == [qQnN] ]]; then
    exit
  fi
done

printf "\n Please make sure that :\nYour data are in the same root folder as this script meaning :\n"
printf "  in the data/ folder there should be :\n"
printf "  - this script \"first_kraken.sh\" \n"
printf "  - kraken_db folder with kraken databases"
printf "  - VANA/trimmed_cutadapt folder with preprocessed VANA data .fastq R1 and R2 paired end files \n"
printf "  - siRNA/trimmed_cutadapt folder with preprocessed siRNA data .fastq files \n"

UserChoice=0
while [[ $UserChoice != [123] ]]
#jsp
do
  echo " ---------------------------------"
	printf "\n\n There are two methods to classify the data and select the potential viral reads.\n"
  printf "  - You can choose to process by elimination, and delete all reads that aren't viral. \n"
  printf "    In this case, you can choose to eliminate every reads matching human, bacteria, archaea and fungi genomes [1]. \n And you can add a plant database to this process, and also delete every read matching plant genomes [2]. \n"
	printf "  - Or you can choose to use kraken viral database and keep the reads that match the database. [3]\n"
	printf "\n Note that all reads, classified and unclassified, will be conserved and accessible.\n \n"
  printf "Which method do you want to use ? \n Please type :\n \"1\" for elimination matching human, bacteria, archaea and fungi genomes\n \"2\" for elimination matching human, bacteria, archae and fungi AND plant genomes\n \"3\" for direct matching on kraken viral database\n \"q\" to quit.\n"
  echo " ---------------------------------"
  read -p ' Method to use : ' UserChoice
	echo " ---------------------------------"
  if [[ $UserChoice == [qQ] ]]; then
    exit
  fi
done

elim1=false
elim2=false
map=false

case  $UserChoice in 

	1)
	  elim1=true
	;;
	2)
		elim1=true
	  elim2=true
	;;
	3)
	  map=true
	;;
	13)
		elim1=true
		map=true
	;;
	12)
		elim1=true
		elim2=true
	;;
	23)
		elim1=true
		elim2=true
		map=true
	;;
	123)
		elim1=true
		elim2=true
		map=true
	;;
esac

printf "$elim1 $elim2 $map"
UserChoice=0
while [[ $UserChoice != [123] ]]
do
  printf " Which data do you want to process ? \n Please type :\n \"1\" for VANA\n \"2\" for siRNA\n \"3\" for both\n \"q\" to quit.\n"
  echo " ---------------------------------"
  read -p ' Data to process : ' UserChoice
  echo " ---------------------------------"
  if [[ $UserChoice == [qQ] ]]; then
    exit
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

############
# VANA

if [ "$VANA" = true ]; then

	# Creation directories :
	cd $SCRIPT_DIR
	echo "Creating VANA directories..."
	trimmed_VANA="$VANA_DIR/trimmed_cutadapt"
	cd $trimmed_VANA

	for fileR1 in $trimmed_VANA/*R1*.fastq
	do
		fileR1=$(basename $fileR1 .fastq)
		cd $trimmed_VANA/$fileR1
		mkdir filtered_reads_viral -p

		if [ "$elim1" = true]; then
			cd $trimmed_VANA/$fileR1/filtered_reads_viral
			mkdir no_viral
			cd no_viral
			#kraken
			kraken2 --threads 10 --db ~/kraken_db/kraken2_no_viral_db --unclassified-out potential_viralR1.fa --classified-out non_viralR1.fa --report report_no_viralR1 --use-names "$trimmed_VANA/$fileR1.fastq"
		fi

		if [ "$elim2" = true]; then
			cd $trimmed_VANA/$fileR1/filtered_reads_viral
			mkdir plant
			cd plant
			#kraken
			kraken2 --threads 10 --db ~/kraken_db/kraken2_plant --unclassified-out potential_viralR1.fa --classified-out plantR1.fa --report report_plantR1 --use-names ../no_viral/potential_viral.fa
		fi

		if [ "$map" = true]; then
			cd $trimmed_VANA/$fileR1/filtered_reads_viral
			mkdir viral
			cd viral
			#kraken
			kraken2 --threads 10 --db ~/kraken_db/kraken2_no_viral_db --unclassified-out potential_viralR1.fa --classified-out non_viralR1.fa --report report_no_viralR1 --use-names "$trimmed_VANA/$fileR1.fastq"
		fi 
	done

	for fileR2 in $trimmed_VANA/*R2*.fastq
	do
		fileR2=$(basename $fileR2 .fastq)
		dirR2=$(basename ${fileR2/R2/R1} .fastq)
		cd $trimmed_VANA/$dirR2
		mkdir filtered_reads_viral -p
		cd filtered_reads_viral

		if [ "$elim1" = true]; then
			cd $trimmed_VANA/$dirR2/filtered_reads_viral
			mkdir no_viral
			cd no_viral
			#kraken
			kraken2 --threads 10 --db ~/kraken_db/kraken2_no_viral_db --unclassified-out potential_viralR2.fa --classified-out non_viralR2.fa --report report_no_viralR2 --use-names "$trimmed_VANA/$fileR2.fastq"
		fi

		if [ "$elim2" = true]; then
			cd $trimmed_VANA/$dirR2/filtered_reads_viral
			mkdir plant
			cd plant
			#kraken
			kraken2 --threads 10 --db ~/kraken_db/kraken2_plant --unclassified-out potential_viralR2.fa --classified-out plantR2.fa --report report_plantR2 --use-names ../no_viral/potential_viral.fa
		fi

		if [ "$map" = true]; then
			cd $trimmed_VANA/$dirR2/filtered_reads_viral
			mkdir viral
			cd viral
			#kraken
			kraken2 --threads 10 --db ~/kraken_db/kraken2_no_viral_db --unclassified-out potential_viralR2.fa --classified-out non_viralR2.fa --report report_no_viralR2 --use-names "$trimmed_VANA/$fileR2.fastq"
		fi 
	done
fi


############
#siRNA

if [ "$siRNA" = true ]; then
	
	# Creation directories :
	cd $SCRIPT_DIR
	echo "Creating siRNA directories..."
	trimmed_siRNA="$siRNA_DIR/trimmed_cutadapt"
	cd $trimmed_siRNA

	for file in $trimmed_siRNA/*.fastq
	do
		fileR1=$(basename $file .fastq)
		cd $trimmed_VANA/$file
		mkdir filtered_reads_viral -p

		if [ "$elim1" = true]; then
			cd $trimmed_VANA/$file/filtered_reads_viral
			mkdir no_viral
			cd no_viral
			#kraken
			kraken2 --threads 10 --db ~/kraken_db/kraken2_no_viral_db --unclassified-out potential_viral.fa --classified-out non_viral.fa --report report_no_viral --use-names "$trimmed_siRNA/fileR1.fastq"
		fi

		if [ "$elim2" = true]; then
			cd $trimmed_VANA/$file/filtered_reads_viral
			mkdir plant
			cd plant
			#kraken
			kraken2 --threads 10 --db ~/kraken_db/kraken2_plant --unclassified-out potential_viral.fa --classified-out plant.fa --report report_plant --use-names ../no_viral/potential_viral.fa
		fi

		if [ "$map" = true]; then
			cd $trimmed_VANA/$file/filtered_reads_viral
			mkdir viral
			cd viral
			#kraken
			kraken2 --threads 10 --db ~/kraken_db/kraken2_no_viral_db --unclassified-out potential_viral.fa --classified-out non_viral.fa --report report_no_viral --use-names "$trimmed_siRNA/$file.fastq"
		fi 
	done

fi

