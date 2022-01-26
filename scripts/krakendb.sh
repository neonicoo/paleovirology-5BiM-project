#!/bin/bash -i

# Authors: Adrian Zurmely - Maelle Broustal - Nicolas Mendiboure INSA Lyon 5BIM

###################### Download,  create and install  kraken databases #######################
##############################################################################################

## donne le repertoire du script : 
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
mkdir kraken_db
KRAKEN_DIR="${SCRIPT_DIR}/kraken_db"

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

# ajouter question etc
printf "This script allows you to install the kraken databases needed to run the pipeline.\n
This operation could take a lot of time."
while [[ $UserChoice != [yn] ]]
do
  printf "Proceed ?"
  read -p '[y/n]' UserChoice
  if [[ $UserChoice == [qQnN] ]]; then
    exit
  fi
done

cd $KRAKEN_DIR
#let's go
echo "Installing NCBI taxonomy..."
kraken2-build --threads 10 --download-taxonomy --db kraken2_no_viral_db &
mkdir kraken2_plant kraken_viral no_viral_18 plant_18 viral18 &
wait
echo "    Taxonomy downloaded"
echo "Copy of taxonomy..."
cp -r kraken2_no_viral_db/taxonomy kraken_viral/taxonomy &
cp -r kraken2_no_viral_db/taxonomy kraken2_plant/taxonomy &
cp -r kraken2_no_viral_db/taxonomy no_viral_18/taxonomy &
cp -r kraken2_no_viral_db/taxonomy plant_18/taxonomy &
cp -r kraken2_no_viral_db/taxonomy viral18/taxonomy &
wait
echo "    Copied."

printf "--------------------------------------------------\n"
echo "Installing library..."
echo "    for VANA no_viral database..."
kraken2-build --threads 10 --download-library  human --db kraken2_no_viral_db &
kraken2-build --threads 10 --download-library  bacteria --db kraken2_no_viral_db &
kraken2-build --threads 10 --download-library  archaea --db kraken2_no_viral_db &
kraken2-build --threads 10 --download-library  fungi --db kraken2_no_viral_db &
wait
echo "       Installed"
#kraken2_plant
echo "    for VANA plant database..."
kraken2-build --threads 10 --download-library plant --db kraken2_plant &
wait
echo "       Installed"
#kraken_viral
echo "    for VANA viral database..."
kraken2-build --threads 10 --download-library viral --db kraken_viral &
wait
echo "       Installed"
#no_viral_18
echo "    for siRNA no_viral database..."
cp -r kraken2_no_viral_db/library no_viral_18/library &
wait
echo "       Installed"
#plant_18
echo "    for siRNA plant database..."
cp -r kraken2_plant/library plant_18/library &
wait
echo "       Installed"
#viral18
echo "    for siRNA viral database..."
cp -r kraken_viral/library viral18/library &
wait
echo "       Installed"

printf "--------------------------------------------------\n"
echo "Building the databases"
echo "May take a long time..."
kraken2-build --threads 10 --build --db kraken2_no_viral_db &
kraken2-build --threads 10 --build --db kraken2_plant &
kraken2-build --threads 10 --build --db kraken_viral &
kraken2-build --threads 10 --build --kmer-len 18 --minimizer-len 11 --minimizer-spaces 2 --db no_viral_18 &
kraken2-build --threads 10 --build --kmer-len 18 --minimizer-len 11 --minimizer-spaces 2 --db  plant_18 &
kraken2-build --threads 10 --build --kmer-len 18 --minimizer-len 11 --minimizer-spaces 2 --db  viral18 &
wait
echo "     Built"
printf "--------------------------------------------------\n"
echo "Cleaning (will help regain storage space) ..."
kraken2-build --threads 10 --clean --db kraken2_no_viral_db &
kraken2-build --threads 10 --clean --db kraken2_plant &
kraken2-build --threads 10 --clean --db kraken_viral &
kraken2-build --threads 10 --clean --db no_viral_18 &
kraken2-build --threads 10 --clean --db  plant_18 &
kraken2-build --threads 10 --clean --db  viral18 &
wait
echo "	  Cleaned"
printf "--------------------------------------------------\n"
conda deactivate
echo "All databases used in this pipeline are constructed and ready to use!"
printf "\n\n\n"
