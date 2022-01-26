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

##################################### Pre-processing #########################################
##############################################################################################


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
	echo "$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
	
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
	echo "$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

	FILES="raw"
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
	for namefileR1 in $FILES/*R1*.fastq
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
		echo "spades.py -o contigs_rnaviral/ -s $FILE -k ${i} --rnaviral "
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

	echo "compress_fasta on "$1
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

printf "\nHello, this is the assembly script for VANA and siRNA data.\n"
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
	FILES_FA="${VANA_DIR}/trimmed_cutadapt"
	mkdir $FILES_FA/SPAdes	-p
	R1="R1"
	R2="R2"
	for namefileR1 in ${FILES_FA}/*R1*.fastq
	do
	  #echo "$(basename $namefileR1 .fastq)"
	  namefileR1=$(basename $namefileR1 .fastq)
	  #on remplace le R1 par le R2 :
	  namefileR2=${namefileR1/R1/R2}
	  #echo "R2 = ${namefileR2/R1/R2}"
	  cd $FILES_FA/
	  echo "$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
	  echo "Creating $namefileR1 directory"
	  mkdir $namefileR1 -p
	  cd $namefileR1
	  echo "Assembling $namefileR1 and $namefileR2 pair-end files..."
	  
	  # on fait le bbmap/repair
	  bash ${BBMAP}/repair.sh \
	  		 in=$FILES_FA/$namefileR1.fastq \
	  		 in2=$FILES_FA/$namefileR2.fastq \
	  		 out1=$FILES_FA/$namefileR1/$namefileR1.repaired_R1.fastq \
	  		 out2=$FILES_FA/$namefileR1/$namefileR2.repaired_R2.fastq \
	  		 outs=$FILES_FA/$namefileR1/$namefileR1.unpaired.fastq
	  
	  #SI $FILES_FA/$namefileR1/$namefileR1.unpaired.fastq est vide : spades va bugger 
	  if [ -s $FILES_FA/$namefileR1/$namefileR1.unpaired.fastq ]
	  then
		# lance spades : 
		spades.py -o outputMETA/ \
							-1 $FILES_FA/$namefileR1/$namefileR1.repaired_R1.fastq \
							-2 $FILES_FA/$namefileR1/$namefileR2.repaired_R2.fastq \ 
							-s $FILES_FA/$namefileR1/$namefileR1.unpaired.fastq \
							--meta
	  else
		spades.py -o outputMETA/ \
							-1 $FILES_FA/$namefileR1/$namefileR1.repaired_R1.fastq \
							-2 $FILES_FA/$namefileR1/$namefileR2.repaired_R2.fastq \
							--meta
	  fi
	  
	  #pour idba_ud on merge 1 et 2 en un seul fichier sans les unpaired :
	  fq2fa --merge $FILES_FA/$namefileR1/$namefileR1.repaired_R1.fastq \
	  							$FILES_FA/$namefileR1/$namefileR2.repaired_R2.fastq \
	  							$FILES_FA/$namefileR1/$namefileR1.mergedR1-R2.fastq

	  idba_ud -r $FILES_FA/$namefileR1/$namefileR1.mergedR1-R2.fastq \
	  				--num_threads 8 \
	  				-o idba_ud_out

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
	FILES_FA="${SIRNA_DIR}/trimmed_cutadapt"
	for f in ${FILES_FA}/*.fastq
	do
	  cd $FILES_FA/
	  f=$(basename $f .fastq)
	  echo "Creating $f directory"
	  mkdir $f -p
	  cd $f
	  echo "$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
	  mkdir SPAdes -p
	  # lance spades_merge_siRNA sur fichier fasta 1e arg, dossier sortie 2e arg, kmers de 13 à 37 pas de 2 :
	  spades_merge $FILES_FA/$f.fastq $FILES_FA/$f/SPAdes/ 13 37 2
	  
	  # compress les fasta avec le script de denis si y a redondance : 
	  cd $FILES_FA/$f/SPAdes/

    cp -p merged_SPAdes.fasta ${f}_compressed_SPAdes.fasta

    compress_myfasta ${f}_compressed_SPAdes.fasta
    mv ${f}_compressed_SPAdes.fasta_compressed.fa ${f}_novel_compressed_SPAdes.fa
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
		  	-evalue 0.001 \
		  	-outfmt 6

	echo "Blastx on virusdetect vrl_plant DB"
	blastx -query $1 \
			-db ${VIRUSDETECTDB_DIR}/vrl_Plants_239_U100_prot\
			-out ${outputx}\
			-num_threads 8 \
			-evalue 0.001 \
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