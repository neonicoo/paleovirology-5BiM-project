#!/bin/bash
#assemblage maison



#--------------------------------------------
# Authors: Adrian Zurmely - Maelle Broustal - Nicolas Mendiboure INSA Lyon 5BIM

#VANA :
#spades.py -o outMETA/ -1 Forward_R1.fastq -2 Reverse_R1.fastq --meta
#siRNA :
#spades.py -o test_kmer/ -s ../../190725_SNK268_B_L006_AFVV-15_R1_trimmed.fastq --rnaviral
#ou sans --rnaviral

# test 
# spades.py -o test_kmer/ -s ../../190725_SNK268_B_L006_AFVV-15_R1_trimmed.fastq -k 11,13,15,17,19,21,23,25,27,29,31,33,35,37 --rnaviral
 

 
# donne le repertoire du script : 
#SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
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



#######
# VANA
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
	

#######
#siRNA
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
	  # lance spades_merge_siRNA sur fichier fasta 1e arg, dossier sortie 2e arg, kmers de 13 à 37 pas de 2 :
	  bash ~/scripts/spades_merge_siRNA.sh $FILES_FA/$f.fastq $FILES_FA/$f/SPAdes/ 13 37 2
	  #bash idba_merge_siRNA.sh . $f 13 37 2
	  
	  # compress les fasta avec le script de denis si y a redondance : 
	  cd $FILES_FA/$f/SPAdes/
      cp merged_SPAdes.fasta compressed_SPAdes.fasta
      bash ~/scripts/compress_fasta.sh  compressed_SPAdes.fasta
	done
fi
	
	
	
	
	
	
	
	
	
	
	
	
