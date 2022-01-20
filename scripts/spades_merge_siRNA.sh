#!/bin/bash
# spades_merge_siRNA.sh


#--------------------------------------------
# Authors: Adrian Zurmely - Maelle Broustal - Nicolas Mendiboure INSA Lyon 5BIM
# ex : spades_merge_siRNA.sh 
# bash ~/scripts/spades_merge_siRNA.sh $FILES_FA/$f.fastq $FILES_FA/$f/SPAdes/ 11 37 2
	  

# donne le dossier en argument 1, le fichier en 2e argument, la taille de kmer minimal en 3e argument, puis la taille de kmer max en 4e argument puis le pas en 5e argument espacés par des espaces.
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







