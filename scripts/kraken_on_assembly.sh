#!/bin/bash
# kraken_on_assembly.sh

#--------------------------------------------
# Authors: Adrian Zurmely - Maelle Broustal - Nicolas Mendiboure INSA Lyon 5BIM

# Une commande ça a ce format là: 
# kraken2 --threads 10 -db kraken_db/viral18 --unclassified-out fichier_unclassified --classified-out fichier_classified --report report --use-names fichier_source
# ça c'est pour tester avec les kmers taille 18
# et pour les kmers taille 15 faut prendre la db kraken_db/viral15, pour les kmers taille 11 c'est kraken_db/little_viral


# Les fichiers en entrée de kraken sur le pipeline principal c'est dans 
# VANA :
# /test_data/VANA/trimmed_cutadapt/nom_echantillon/outputMETA/contigs.fasta
# /test_data/VANA/trimmed_cutadapt/nom_echantillon/idba_ud_out/contig.fa
# et pour siRNA c'est un peu plus technique c'est dans :
# /test_data/siRNA/trimmed_cutadapt/nom_echantillon/SPAdes/compressed_SPAdes.fasta_Compressed.fa


DOSSIER_siRNA=~/test_data/siRNA/trimmed_cutadapt/
DOSSIER_VANA=~/data/VANA/trimmed_cutadapt/


for namefileR1 in $DOSSIER_VANA/*R1*.fastq
do 
	namefileR1=$(basename $namefileR1 .fastq)
	cd $DOSSIER_VANA/$namefileR1
	mkdir kraken -p
	cd kraken
	
	mkdir SPAdes -p
	cd SPAdes
	mkdir 1_no_viraldb
	mkdir 2_plant
	mkdir 3_viral
# kraken sur recup outputMETA/contigs.fasta
	# 1e etape
	kraken2 --threads 8 -db ~/kraken_db/kraken2_no_viral_db --unclassified-out 1_no_viraldb/potential_viral.fa --classified-out 1_no_viraldb/non_viral.fa --report 1_no_viraldb/report_no_viral --use-names $DOSSIER_VANA/$namefileR1/outputMETA/contigs.fasta
	# 2e etape
	kraken2 --threads 8 -db ~/kraken_db/kraken2_plant --unclassified-out 2_plant/potential_viral.fa --classified-out 2_plant/plant.fa --report 2_plant/report_plant --use-names $DOSSIER_VANA/$namefileR1/kraken/SPAdes/1_no_viraldb/potential_viral.fa
	# 3e etape 
	kraken2 --threads 8 -db ~/kraken_db/kraken_viral --unclassified-out 3_viral/non_viral.fa --classified-out 3_viral/viral.fa --report 3_viral/report_viral --use-names $DOSSIER_VANA/$namefileR1/outputMETA/contigs.fasta
	
	mkdir IDBA -p
	cd IDBA
	mkdir 1_no_viraldb
	mkdir 2_plant
	mkdir 3_viral	
# kraken sur recup idba_ud_out/contig.fa
	# 1e etape
	kraken2 --threads 8 -db ~/kraken_db/kraken2_no_viral_db --unclassified-out 1_no_viraldb/potential_viral.fa --classified-out 1_no_viraldb/non_viral.fa --report 1_no_viraldb/report_no_viral --use-names $DOSSIER_VANA/$namefileR1/idba_ud_out/contig.fa
	# 2e etape
	kraken2 --threads 8 -db ~/kraken_db/kraken2_plant --unclassified-out 2_plant/potential_viral.fa --classified-out 2_plant/plant.fa --report 2_plant/report_plant --use-names $DOSSIER_VANA/$namefileR1/kraken/IDBA/1_no_viraldb/potential_viral.fa
	# 3e etape 
	kraken2 --threads 8 -db ~/kraken_db/kraken_viral --unclassified-out 3_viral/non_viral.fa --classified-out 3_viral/viral.fa --report 3_viral/report_viral --use-names $DOSSIER_VANA/$namefileR1/idba_ud_out/contig.fa
	
done



for namefileR1 in $DOSSIER_siRNA/*.fastq
do 
	namefileR1=$(basename $namefileR1 .fastq)
	cd $DOSSIER_siRNA/$namefileR1
	mkdir kraken -p
	cd kraken
	
	mkdir SPAdes -p
	cd SPAdes
	mkdir 1_no_viraldb
	mkdir 2_plant
	mkdir 3_viral
# kraken sur recup 
	# 1e etape
	kraken2 --threads 8 -db ~/kraken_db/no_viral_18 --unclassified-out 1_no_viraldb/potential_viral.fa --classified-out 1_no_viraldb/non_viral.fa --report 1_no_viraldb/report_no_viral --use-names $DOSSIER_siRNA/$namefileR1/SPAdes/compressed_SPAdes.fasta_Compressed.fa
	# 2e etape
	kraken2 --threads 8 -db ~/kraken_db/plant_18 --unclassified-out 2_plant/potential_viral.fa --classified-out 2_plant/plant.fa --report 2_plant/report_plant --use-names $DOSSIER_siRNA/$namefileR1/kraken/SPAdes/1_no_viraldb/potential_viral.fa
	# 3e etape 
	kraken2 --threads 8 -db ~/kraken_db/viral18 --unclassified-out 3_viral/non_viral.fa --classified-out 3_viral/viral.fa --report 3_viral/report_viral --use-names $DOSSIER_siRNA/$namefileR1/SPAdes/compressed_SPAdes.fasta_Compressed.fa
	
done




