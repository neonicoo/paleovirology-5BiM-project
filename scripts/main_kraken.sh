#!/bin/bash -i

# Authors: Adrian Zurmely - Maelle Broustal - Nicolas Mendiboure INSA Lyon 5BIM

##################################### Initial settings  ######################################
##############################################################################################

## donne le repertoire du script : 
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
#echo "$SCRIPT_DIR"
DATA_DIR="${SCRIPT_DIR}/data"
SIRNA_DIR="$DATA_DIR/siRNA"
VANA_DIR="$DATA_DIR/VANA"
DATABASES_DIR="${SCRIPT_DIR}/databases"
BBMAP="${SCRIPT_DIR}/bbmap"
KRAKEN_DIR="${SCRIPT_DIR}/kraken_db"
BLAST_DIR="${DATA_DIR}/blast"
SRC_DIR="${SCRIPT_DIR}/src"
VIRUSDETECTDB_DIR="${DATABASES_DIR}/plant_239_U100"
NR_DIR="${DATABASES_DIR}/nr"
NT_DIR="${DATABASES_DIR}/nt"


if [ -d ${BBMAP} ]
then
	:
else
	tar -xzvf bbmap.tar.gz
	rm -rf bbmap.tar.gz
fi

if [ -d ${DATABASES_DIR} ]
then
	:
else
	mkdir databases 
fi

cd databases

if [ -d $VIRUSDETECTDB_DIR ] 
then
    :
else
    mkdir plant_239_U100
	cd plant_239_U100

	wget bioinfo.bti.cornell.edu/ftp/program/VirusDetect/virus_database/v239/U100/plant_239_U100.tar.gz
	tar -xzvf plant_239_U100.tar.gz --strip-components=1

	wget bioinfo.bti.cornell.edu/ftp/program/VirusDetect/virus_database/v239/vrl_genbank.info.gz
	gunzip vrl_genbank.info.gz 
	wget bioinfo.bti.cornell.edu/ftp/program/VirusDetect/virus_database/v239/vrl_idmapping.gz
	gunzip vrl_idmapping.gz
	cd ..
fi

if [ -d $NR_DIR ] 
then
    :
else
    printf "NT database is missing, would you like to download it ? \n"
    printf "This might take quite a few time\n"
    
    UserChoice=0
	while [[ $UserChoice != [12] ]]
	do
	  echo "---------------------------------"
	  printf "Please type : \n \"1\" for YES\n \"2\" for NO\n \"q\" to quit.\n"
	  read -p 'Dowload NR DB for blast ? : ' UserChoice
	  if [[ $UserChoice == [qQ] ]]; then
	    break
	  fi
	done
	if [[ $UserChoice == 1 ]]
	then
		mkdir $NR_DIR
		cd $NR_DIR
		update_blastdb nr
	else
		break
	fi
fi

if [ -d $NT_DIR ] 
then
    :
else
    printf "NT database is missing, would you like to download it ? \n"
    printf "This might take quite a few time\n"
    
    UserChoice=0
	while [[ $UserChoice != [12] ]]
	do
	  echo "---------------------------------"
	  printf "Please type : \n \"1\" for YES\n \"2\" for NO\n \"q\" to quit.\n"
	  read -p 'Dowload NT DB for blast ? : ' UserChoice
	  if [[ $UserChoice == [qQ] ]]; then
	    break
	  fi
	done
	if [[ $UserChoice == 1 ]]
	then
		mkdir $NT_DIR
		cd $NT_DIR
		update_blastdb nt
	else
		break
	fi
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
 

if [ -d ${SIRNA_DIR} ] 
then
    :
else
    echo "${SIRNA_DIR} doesn't exist"
    exit
fi

if [ -d ${VANA_DIR} ] 
then
    :
else
    echo "${VANA_DIR} doesn't exist"
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
	cd $SCRIPT_DIR
	echo "Preprocessing siRNA data"
	cd $SIRNA_DIR
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
	cd $VANA_DIR
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
        cd $VANA_DIR

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


###################################### Filtrage kraken ######################################
##############################################################################################

if [ -d "${SIRNA_DIR}/trimmed_cutadapt" ] 
then
    :
else
    echo "${SIRNA_DIR}/trimmed_cutadapt doesn't exist"
    break
fi

if [ -d "${VANA_DIR}/trimmed_cutadapt" ] 
then
    :
else
    echo "${VANA_DIR}/trimmed_cutadapt doesn't exist"
    break
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
    break
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
	printf "\n There are two methods to classify the data and select the potential viral reads.\n"
  printf "  - You can choose to process by elimination, and delete all reads that aren't viral. \n"
  printf "    In this case, you can choose to eliminate every reads matching human, bacteria, archaea and fungi genomes [1]. \n And you can add a plant database to this process, and also delete every read matching plant genomes [2]. \n"
	printf "  - Or you can choose to use kraken viral database and keep the reads that match the database. [3]\n"
	printf "\n Note that all reads, classified and unclassified, will be conserved and accessible.\n \n"
  printf "Which method do you want to use ? \n Please type :\n \"1\" for elimination matching human, bacteria, archaea and fungi genomes\n \"2\" for elimination matching human, bacteria, archae and fungi AND plant genomes\n \"3\" for direct matching on kraken viral database\n \"q\" to quit.\n"
  echo " ---------------------------------"
  read -p ' Method to use : ' UserChoice
	echo " ---------------------------------"
  if [[ $UserChoice == [qQ] ]]; then
    break
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

UserChoice=0
while [[ $UserChoice != [123] ]]
do
  printf " Which data do you want to process ? \n Please type :\n \"1\" for VANA\n \"2\" for siRNA\n \"3\" for both\n \"q\" to quit.\n"
  echo " ---------------------------------"
  read -p ' Data to process : ' UserChoice
  echo " ---------------------------------"
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
		file=$(basename $fileR1 .fastq)
		mkdir ${file}
		cd $file
		echo "$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
		mkdir filtered_reads_viral -p
		cd filtered_reads_viral

		if [[ "$elim1" = true ]]; then
			echo "$file: Create no_viral dir"
			mkdir no_viral
			cd no_viral
			#kraken
			kraken2 --threads 10 --db ${KRAKEN_DIR}/kraken2_no_viral_db --unclassified-out potential_viralR1.fa --classified-out non_viralR1.fa --report report_no_viralR1 --use-names "$trimmed_VANA/$fileR1.fastq"
			cd ..
		fi

		if [ "$elim2" = true]; then
			echo "$file: Create plant dir"
			mkdir plant
			cd plant
			#kraken
			kraken2 --threads 10 --db ${KRAKEN_DIR}/kraken2_no_viral_db --unclassified-out potential_viralR1.fa --classified-out non_viralR1.fa --report report_no_viralR1 --use-names ../no_viral/potential_viral.fa
			cd ..
		fi

		if [ "$map" = true]; then
			echo "$file: Create viral dir"
			mkdir viral
			cd viral
			#kraken
			kraken2 --threads 10 --db ${KRAKEN_DIR}/kraken2_no_viral_db --unclassified-out potential_viralR1.fa --classified-out non_viralR1.fa --report report_no_viralR1 --use-names "$trimmed_VANA/$fileR1.fastq"
			cd ..
		fi 
	done

	for fileR2 in $trimmed_VANA/*R2*.fastq
	do
		file=$(basename $fileR2 .fastq)
		dirR2=$(basename ${file/R2/R1} .fastq)
		cd $trimmed_VANA/$dirR2
		mkdir filtered_reads_viral -p
		cd filtered_reads_viral

		if [ "$elim1" = true]; then
			echo "$file: Create no_viral dir"
			mkdir no_viral
			cd no_viral
			#kraken
			kraken2 --threads 10 --db ${KRAKEN_DIR}//kraken2_no_viral_db --unclassified-out potential_viralR2.fa --classified-out non_viralR2.fa --report report_no_viralR2 --use-names "$trimmed_VANA/$file.fastq"
			cd ..
		fi

		if [ "$elim2" = true]; then
			echo "$file: Create plant dir"
			mkdir plant
			cd plant
			#kraken
			kraken2 --threads 10 --db ${KRAKEN_DIR}//kraken2_plant --unclassified-out potential_viralR2.fa --classified-out plantR2.fa --report report_plantR2 --use-names ../no_viral/potential_viral.fa
			cd ..
		fi

		if [ "$map" = true]; then
			echo "$file: Create viral dir"
			mkdir viral
			cd viral
			#kraken
			kraken2 --threads 10 --db ${KRAKEN_DIR}/kraken2_no_viral_db --unclassified-out potential_viralR2.fa --classified-out non_viralR2.fa --report report_no_viralR2 --use-names "$trimmed_VANA/$file.fastq"
			cd ..
		fi 
	done
fi

############
#siRNA

if [ "$siRNA" = true ]; then
	
	# Creation directories :
	cd $SCRIPT_DIR
	echo "Creating siRNA directories..."
	trimmed_siRNA="$SIRNA_DIR/trimmed_cutadapt"
	cd $trimmed_siRNA

	for fileR1 in $trimmed_siRNA/*.fastq
	do
		file=$(basename $fileR1 .fastq)
		mkdir ${file}
		cd $file
		echo "$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
		mkdir filtered_reads_viral -p
		cd filtered_reads_viral

		if [[ "$elim1" = true ]]
		then
			echo "Create no_viral dir"
			mkdir no_viral
			cd no_viral
			#kraken
			kraken2 --threads 10 --db ${KRAKEN_DIR}/kraken2_no_viral_db --unclassified-out potential_viral.fa --classified-out non_viral.fa --report report_no_viral --use-names "$trimmed_siRNA/fileR1.fastq"
			cd ..
		fi

		if [[ "$elim2" = true ]]
		then
			echo "Create plant_dir"
			mkdir plant
			cd plant
			#kraken
			kraken2 --threads 10 --db ${KRAKEN_DIR}/kraken2_plant --unclassified-out potential_viral.fa --classified-out plant.fa --report report_plant --use-names ../no_viral/potential_viral.fa
			cd ..
		fi

		if [[ "$map" = true ]]
		then
			echo "Create viral dir"
			mkdir viral
			cd viral
			#kraken
			kraken2 --threads 10 --db ${KRAKEN_DIR}/kraken_viral_db --unclassified-out potential_viral.fa --classified-out non_viral.fa --report report_no_viral --use-names "$trimmed_siRNA/$file.fastq"
			cd ..
		fi 
	done

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
	# Usage: ./compress_fasta.sh File.fa
	# Remove duplicated sequences in a multifasta file
	# Be careful: sequence IDs must be different
	# Output: $1_compressed.fa

	module load fastx_toolkit/0.0.14
	module load vsearch/2.14.0

	echo "./compress_fasta.sh "$1
	echo $(grep -c "^>" $1)" sequences in "$1

	vsearch --cluster_fast $1 --centroids $1_compressed.fa --iddef 0 --id 1.00 --strand both --qmask none --fasta_width 0 --minseqlength 1 --maxaccept 0 --maxreject 0

	fasta_formatter -w 0 -i $1 -o $1.tab -t
	fasta_formatter -w 0 -i $1_compressed.fa -o $1_compressed.fa.tab -t
	cut -f2 $1_compressed.fa.tab | rev | tr atgcATGC tacgTACG > $1_compressed.fa.tab.tab
	paste $1_compressed.fa.tab $1_compressed.fa.tab.tab > $1_compressed.fa.tab.tab.tab
	rm $1_compressed.fa.tab
	rm $1_compressed.fa.tab.tab
	awk -F $"\t" '{print ">"$1"|RC\n"$3}' $1_compressed.fa.tab.tab.tab > $1_compressed.fa.rc 
	rm $1_compressed.fa.tab.tab.tab
	cat $1_compressed.fa $1_compressed.fa.rc > $1_compressed.fa.all
	rm $1_compressed.fa.rc
	old_IFS=$IFS
	IFS=$'\t'
	> $1_compressed.fa.tab
	while read c1 c2
		do
		printf "$c1\t$c2\t" >> $1_compressed.fa.tab
		grep -c $c2 $1_compressed.fa.all >> $1_compressed.fa.tab
	done < $1.tab
	IFS=$old_IFS
	rm $1.tab
	rm $1_compressed.fa.all
	awk -F $"\t" '{if ($3==0) print ">"$1"\n"$2}' $1_compressed.fa.tab > $1_compressed.fa.more
	rm $1_compressed.fa.tab
	echo $(grep -c "^>" $1_compressed.fa.more)" readded sequences"
	cat $1_compressed.fa.more >> $1_compressed.fa
	rm $1_compressed.fa.more 

	fasta_formatter -w 0 -i $1_compressed.fa -o $1_compressed.fa.tab -t 
	awk -F $"\t" '{print $1"\t"$2"\t"length($2)}' $1_compressed.fa.tab | sort -t $'\t' -n -r -k3,3 | awk -F $"\t" '{print ">"$1"\n"$2}' > $1_compressed.fa
	rm $1_compressed.fa.tab

	echo $(grep -c "^>" $1_compressed.fa)" sequences in "$1_compressed.fa
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

printf "\nHello, this is the assembly script for VANA and siRNA data.\n\n Please make sure that you have previously obtained preprocessed data. \n"


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
#Choice between preprocessed and preprocessed-krakened

if [ "$VANA" = true ]; then
	
	#Vérifie que les données filtrées par kraken existent
	krakened=true
	for file in "$VANA_DIR/trimmed_cutadapt"
	do
		if [[ -d "${VANA_DIR}/trimmed_cutadapt/$file/filtered_reads_viral" ]]
		then
			:
		else 
			krakened=false
		fi
	done

	#Choix des données à utiliser: si krakened, choix; sinon, processed anyway
	if [ "$krakened" = true ]; then
		UserChoice=0
		while [[ $UserChoice != [12] ]]
		do
		  echo "---------------------------------"
		  printf "We see that you have filter the preprocess VANA data using kraken.\n You can now choose to assemble the preprocessed data, or the preprocessed and filtered data. \n Which data do you want to assemble ? \n Please type :\n \"1\" for preprocessed data \n \"2\" for prepocessed and filtered data\n \"3\" for both\n \"q\" to quit.\n"
		  read -p 'Data to preprocess : ' UserChoice
		  if [[ $UserChoice == [qQ] ]]; then
			break
		  fi
		done
		
		processed=false
		filtered=false
		case  $UserChoice in 

			3)
			  processed=true
			  filtered=true
			;;
			2)
			  filtered=true
			;;
			1)
			  processed=true
			;;
		esac
	else
		processed=true
	fi
	
	if [ "$processed" = true ]; then

		echo "Assembling preprocessed VANA data"
		FILES_FA="${VANA_DIR}/trimmed_cutadapt"
		mkdir $FILES_FA/outputMETA	-p
		R1="R1"
		R2="R2"
		for namefileR1 in $FILES_FA/*R1*.fastq
		do
			#echo "$(basename $namefileR1 .fastq)"
			namefileR1=$(basename $namefileR1 .fastq)
			#on remplace le R1 par le R2 :
			namefileR2=${namefileR1/R1/R2}
			#echo "R2 = ${namefileR2/R1/R2}"
			cd $FILES_FA/$namefileR1
			echo "Assembling $namefileR1 and $namefileR2 pair-end files..."
			  
			# on fait le bbmap/repair
			bash ${BBMAP}/repair.sh in=$FILES_FA/$namefileR1.fastq in2=$FILES_FA/$namefileR2.fastq out1=$FILES_FA/$namefileR1/$namefileR1.repaired_R1.fastq out2=$FILES_FA/$namefileR1/$namefileR2.repaired_R2.fastq outs=$FILES_FA/$namefileR1/$namefileR1.unpaired.fastq
		  
			#SI $FILES_FA/$namefileR1/$namefileR1.unpaired.fastq est vide : spades va bugger 
			if [ -s $FILES_FA/$namefileR1/$namefileR1.unpaired.fastq ]; then
				# lance spades : 
				spades.py -o outputMETA/ -1 $FILES_FA/$namefileR1/$namefileR1.repaired_R1.fastq -2 $FILES_FA/$namefileR1/$namefileR2.repaired_R2.fastq -s $FILES_FA/$namefileR1/$namefileR1.unpaired.fastq --meta
		  	else
				spades.py -o outputMETA/ -1 $FILES_FA/$namefileR1/$namefileR1.repaired_R1.fastq -2 $FILES_FA/$namefileR1/$namefileR2.repaired_R2.fastq --meta
		  	fi
		  
		  	#pour idba_ud on merge 1 et 2 en un seul fichier sans les unpaired :
		  	fq2fa --merge $FILES_FA/$namefileR1/$namefileR1.repaired_R1.fastq $FILES_FA/$namefileR1/$namefileR2.repaired_R2.fastq $FILES_FA/$namefileR1/$namefileR1.mergedR1-R2.fastq
		  	idba_ud -r $FILES_FA/$namefileR1/$namefileR1.mergedR1-R2.fastq --num_threads 8 -o idba_ud_out

	  	done
	fi
	
	if [ "$filtered" = true ]; then

		echo "Assembling preprocessed VANA data"
		FILES_BASE="${VANA_DIR}/trimmed_cutadapt"
		#mkdir $FILES_FA/outputMETA	-p
		R1="R1"
		R2="R2"
		for namefileR1 in $FILES_BASE/*R1*.fastq
		do
			dirR1=$(basename $namefileR1)
			if [[ -d "${FILES_BASE}/$dirR1/filtered_reads_viral/" ]]; then
				for meth in ${FILES_BASE}/$dirR1/filtered_reads_viral/
				do
					cd ${FILES_BASE}/$dirR1/filtered_reads_viral/$meth
					 
					mkdir SPAdes_filtered
					namefileR1="potential_viralR1"
					namefileR2="potential_viralR2"
					namefile="potential_viral"
					echo "Assembling $namefileR1 and $namefileR2 pair-end files..."
					  
					# on fait le bbmap/repair
					bash ${BBMAP}/repair.sh in=$namefileR1.fastq in2=$namefileR2.fastq out1=$namefile.repaired_R1.fastq out2=$namefile.repaired_R2.fastq outs=$namefile.unpaired.fastq
				  
					#SI $FILES_FA/$namefileR1/$namefileR1.unpaired.fastq est vide : spades va bugger 
					if [ -s SPAdes_filtered/$namefile.unpaired.fastq ]; then
						# lance spades : 
						spades.py -o SPAdes_filtered/ -1 $namefile.repaired_R1.fastq -2 $namefile.repaired_R2.fastq -s $namefile.unpaired.fastq --meta
				  	else
						spades.py -o SPAdes_filtered/ -1 $namefile.repaired_R1.fastq -2 $namefile.repaired_R2.fastq --meta
				  	fi
				  
				  	#pour idba_ud on merge 1 et 2 en un seul fichier sans les unpaired :
				  	fq2fa --merge/$namefile.repaired_R1.fastq $namefile.repaired_R2.fastq $namefile.mergedR1-R2.fastq 
				  	idba_ud -r $namefile.mergedR1-R2.fastq --num_threads 8 -o idba_ud_out
				done
			fi
	  	done
	fi
fi

#####################
# siRNA #############
#####################

#pour tous les échantillons  siRNA/trimmed_cutadapt/
	# pour tous les fichiers test_data/siRNA/trimmed_cutadapt/SPAdes/BO1_15/test_kmer/K*/final_contigs.fasta :
	# echo le fichier >> merged_fasta_kmer.fasta
	
	
if [ "$siRNA" = true ]; then
	
	#Vérifie que les données filtrées par kraken existent
	krakened=true
	for file in "$siRNA_DIR/trimmed_cutadapt"
	do
		if [[ -d "${siRNA_DIR}/trimmed_cutadapt/$file/filtered_reads_viral" ]]
		then
			:
		else 
			krakened=false
		fi
	done

	#Choix des données à utiliser: si krakened, choix; sinon, processed anyway
	if [ "$krakened" == true ]; then
		UserChoice=0
		while [[ $UserChoice != [12] ]]
		do
		  echo "---------------------------------"
		  printf "We see that you have filter the preprocess siRNA data using kraken.\n You can now choose to assemble the preprocessed data, or the preprocessed and filtered data. \n Which data do you want to assemble ? \n Please type :\n \"1\" for preprocessed data \n \"2\" for prepocessed and filtered data\n \"3\" for both\n \"q\" to quit.\n"
		  read -p 'Data to preprocess : ' UserChoice
		  if [[ $UserChoice == [qQ] ]]; then
			break
		  fi
		done
		
		processed=false
		filtered=false
		case  $UserChoice in 

			3)
			  processed=true
			  filtered=true
			;;
			2)
			  filtered=true
			;;
			1)
			  processed=true
			;;
		esac
	else
		processed=true
	fi
	
	if [ "$processed" = true ]; then
		
		echo "Assembling siRNA data"
		FILES_FA="${SIRNA_DIR}/trimmed_cutadapt/"
		for f in $FILES_FA*.fastq
		do
		  cd $FILES_FA/
		  f=$(basename $f .fastq)
		  cd $f
		  mkdir SPAdes -p
		  # lance spades_merge_siRNA sur fichier fasta 1e arg, dossier sortie 2e arg, kmers de 13 à 37 pas de 2 :
		  spades_merge $FILES_FA/$f.fastq $FILES_FA/$f/SPAdes/ 13 37 2
		  #bash idba_merge_siRNA.sh . $f 13 37 2
		  
		  # compress les fasta avec le script de denis si y a redondance : 
		  cd $FILES_FA/$f/SPAdes/
		  cp merged_SPAdes.fasta compressed_SPAdes.fasta
		  compress_myfasta compressed_SPAdes.fasta
		  mv compressed_SPAdes.fasta_Compressed.fa compressed_SPAdes.fa
		done
	fi
	
	if [ "$filtered" = true ]; then
		
		echo "Assembling siRNA data"
		FILES_BASE="${SIRNA_DIR}/trimmed_cutadapt/"
		for f in $FILES_BASE*.fastq
		do
			cd $FILES_BASE/
			f=$(basename $f .fastq)
			cd $f
			if [[ -d "${FILES_BASE}/$f/filtered_reads_viral/" ]]; then
				for meth in ${FILES_BASE}/$f/filtered_reads_viral/
				do
					cd ${FILES_BASE}/$f/filtered_reads_viral/$meth
					mkdir SPAdes -p
					# lance spades_merge_siRNA sur fichier fasta 1e arg, dossier sortie 2e arg, kmers de 13 à 37 pas de 2 :
					spades_merge $FILES_FA/$f/filtered_reads_viral/$meth/potential_viral.fa $FILES_FA/$f/filtered_reads_viral/$meth/SPAdes/ 13 37 2
					#bash idba_merge_siRNA.sh . $f 13 37 2
					  
					# compress les fasta avec le script de denis si y a redondance : 
					cd $FILES_FA/$f/filtered_reads_viral/$meth/SPAdes/
					cp merged_SPAdes.fasta compressed_SPAdes.fasta
					compress_myfasta compressed_SPAdes.fasta
					mv compressed_SPAdes.fasta_Compressed.fa compressed_SPAdes.fa	
			fi
		done
	fi

########################################## Blast #############################################
##############################################################################################

if [ -d ${BLAST_DIR} ]
then
	:
else
	mkdir ${BLAST_DIR}
	mkdir ${BLAST_DIR}/siRNA
	mkdir ${BLAST_DIR}/VANA
	echo "Created $BLAST_DIR"
fi



function virusdetect_blast ()
{
	echo "$1"

	outputn=${1/_trimmed/}_virusdetect_blastn.txt
	outputx=${1/_trimmed/}_virusdetect_blastx.txt

	echo "Blastn on virusdetect vrl_plant DB"
	blastn -query $1 \
		  	-db ${VIRUSDETECTDB_DIR}/vrl_Plants_239_U100 \
		  	-out ${outputn}\
		  	-num_threads 8 \
		  	-evalue 0.01 \
		  	-outfmt 6

	echo "Blastx on virusdetect vrl_plant DB"
	blastx -query $1 \
			-db ${VIRUSDETECTDB_DIR}/vrl_Plants_239_U100_prot\
			-out ${outputx}\
			-num_threads 8 \
			-evalue 0.01 \
			-outfmt 6

	
	if [[ -s ${outputn} ]]
	then 
		echo "Virus identification from blastn $1 with python3 scripts"
		python3 ${SRC_DIR}/blastn_virus_identify.py \
						${outputn} \
						$VIRUSDETECTDB_DIR/ \
						${outputn/_virusdetect_blastn.txt/}_virusdetect_blastn_taxon
	fi
	if [[ -s ${outputx} ]]
	then 
		echo "Virus identification from blastx $1 with python3 scripts"
		python3 ${SRC_DIR}/blastx_virus_identify.py \
						${outputx} \
						$VIRUSDETECTDB_DIR/ \
						${outputx/_virusdetect_blastx.txt/}_virusdetect_blastx_taxon
	fi
}

function nrnt_blast ()
{

	outputn=${1/_trimmed/}_nt_result_blastn.txt
	outputx=${1/_trimmed/}_nr_result_blastx.txt

	echo "Blastx on nr DB"
	blastx -query $1 \
			-db ${DATABASES_DIR}/nr \
			-out $outputx \
			-num_threads 12 \
			-evalue 0.0001 \
			-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames sskingdoms stitle" \
			-max_target_seqs 2

	echo "Blastn on nt DB"
	blastn -query $1 \
			-db ${DATABASES_DIR}/nt \
			-out $outputn \
			-num_threads 12 \
			-evalue 0.0001 \
			-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames sskingdoms stitle" \
			-max_target_seqs 2
			
	# je fais un grep pour récupérer les blast que de "Viruses"
  	grep -i "Viruses" $outputn.txt  > viruses_BlastN.txt
	# je récupère les identifiants des contigs de virus
  	cut -f1 viruses_BlastN.txt > virus_ID_BlastN.txt
	# on ne garde dans le .fasta que les contigs identifiés comme virus : 
  	grep -f virus_ID_BlastN.txt -i $1 -A 1 > viruses_contigs_BlastN.fasta
	
}


printf "\nHello, this is the blast step.\n"
UserChoice=0
while [[ $UserChoice != [123] ]]
do
  echo "---------------------------------"
  printf "Which data do you want to blast? \n Please type :\n \"1\" for VANA\n \"2\" for siRNA\n \"3\" for both\n \"q\" to quit.\n"
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

if [ "$siRNA" = true ]
then
	for sample_dir in ${SIRNA_DIR}/trimmed_cutadapt/*_trimmed
	do
		cd $sample_dir
		#echo "$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
		sample=$(basename $sample_dir)
		echo "Copy the SPAdes assembled siRNA $sample"
		cd $sample_dir/SPAdes/ 
		cp  -p ${sample}_novel_compressed_SPAdes.fa ${BLAST_DIR}/siRNA/
	done

	UserChoice=0
	while [[ $UserChoice != [123] ]]
	do 
		echo "---------------------------------"
		printf "Which database you want to give to blast your siRNA contigs ? \n Please type :\n \"1\" for Virusdetect vrl_plant DB\n \"2\" for nr and nt DB\n \"3\" for both nrnt and virusdetect \n "
	  	read -p 'Which reference database  : ' UserChoice
	  	if [[ $UserChoice == [qQ] ]]
	  	then
	    	break
	  	fi
	done

	virusdetect=false
	nrnt=false

	case $UserChoice in 
		1)
		 	virusdetect=true
		;;
		2)
		 	nrnt=true
		;;
		3)
			virusdetect=true
			nrnt=true
		;;
	esac

	if [ "$virusdetect" = true ]
	then
		for file in ${BLAST_DIR}/siRNA/*.fa
		do
			virusdetect_blast $file
		done

		mkdir ${BLAST_DIR}/siRNA/virusdetect_results
		mv ${BLAST_DIR}/siRNA/*_virusdetect_blast{n,x}_taxon.txt ${BLAST_DIR}/siRNA/virusdetect_results
		cat ${BLAST_DIR}/siRNA/virusdetect_results/*taxon.txt > vrl_plant_results_siRNA_all.txt
	fi

	if [ "$nrnt" = true ]
	then
		for file in ${BLAST_DIR}/siRNA/*.fa
		do
			nrnt_blast $file
		done
	fi

fi


#####################
# VANA ##############
#####################

if [ "$VANA" = true ]
then
	for sample_dir in ${VANA_DIR}/trimmed_cutadapt/*_trimmed
	do
		cd $sample_dir
		echo "$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
		sample=$(basename $sample_dir)
		echo "Copy the IDBA assembled VANA $sample"
		cd $sample_dir/idba_ud_out/ 
		cp  -p contig.fa ${BLAST_DIR}/VANA/${sample}_idba_ud_contigs.fa
		echo "Copy the SPAdes META assembled VANA $sample"
		cd $sample_dir/outputMETA/
		cp  -p contigs.fasta ${BLAST_DIR}/VANA/${sample}_spades_meta_contigs.fasta
	done

	UserChoice=0
	while [[ $UserChoice != [123] ]]
	do 
		echo "---------------------------------"
		printf "Which database you want to give to blast your VANA contigs ? \n Please type :\n \"1\" for Virusdetect vrl_plant DB\n \"2\" for nr and nt DB\n \"3\" for both nrnt and virusdetect \n "
	  	read -p 'Which reference database  : ' UserChoice
	  	if [[ $UserChoice == [qQ] ]]
	  	then
	    	break
	  	fi
	done

	virusdetect=false
	nrnt=false

	case $UserChoice in 
		1)
		 	virusdetect=true
		;;
		2)
		 	nrnt=true
		;;
		3)
			virusdetect=true
			nrnt=true
		;;
	esac

	if [ "$virusdetect" = true ]
	then
		for file in ${BLAST_DIR}/VANA/*.fa*
		do
			virusdetect_blast $file
		done

		mkdir ${BLAST_DIR}/VANA/virusdetect_results
		mv ${BLAST_DIR}/VANA/*_virusdetect_blast{n,x}_taxon.txt ${BLAST_DIR}/VANA/virusdetect_results
		cat ${BLAST_DIR}/VANA/virusdetect_results/*taxon.txt > vrl_plant_results_VANA_all.txt
	fi

	if [ "$nrnt" = true ]
	then
		for file in ${BLAST_DIR}/VANA/*.fa
		do
			nrnt_blast $file
		done
	fi
fi

conda deactivate
echo "paleogenomic conda env OFF"
