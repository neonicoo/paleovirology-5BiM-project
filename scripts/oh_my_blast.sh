#!/usr/bin/bash -i

FILEPATH=$2 #contigs to blast
FILE="${FILEPATH##*/}"
FILENAME="${FILE%.*}"

conda activate paleogenomic
echo "paleogenomic conda env ON"

while [ -n "$1" ]; do # while loop starts

	case "$1" in

	-virusdetect) 
		echo "You're going to blast{n,x} your contigs over the virus plants database"
		echo -n "Continue? [y/n] " 
		read doit
			case $doit in 
				y|Y) 
					# $3 : virusdetect_db

					blastn -query $FILEPATH \
						  	-db $3vrl_Plants_239_U100 \
						  	-out ${FILENAME}_virusdetect_blastn.txt\
						  	-num_threads 8 \
						  	-evalue 0.01 \
						  	-outfmt 6

					echo "#### blastN on virusdetect db - success ####"

					blastx -query $FILEPATH \
							-db $3vrl_Plants_239_U100_prot\
							-out ${FILENAME}_virusdetect_blastx.txt\
							-num_threads 8 \
							-evalue 0.01 \
							-outfmt 6

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

					blastn -query $FILEPATH \
						  	-db $3 \
						  	-out ${FILENAME}_nr_blastn.txt\
						  	-num_threads 8 \
						  	-evalue 0.01 \
						  	-outfmt 6

					echo "#### blastN on nr db - success ####"

					blastx -query $FILEPATH \
							-db $3 \
							-out ${FILENAME}_nr_blastx.txt\
							-num_threads 8 \
							-evalue 0.01 \
							-outfmt 6
					echo "#### blastX nr db - success ####"

					blastn -query $FILEPATH \
						  	-db $4 \
						  	-out ${FILENAME}_nt_blastn.txt\
						  	-num_threads 8 \
						  	-evalue 0.01 \
						  	-outfmt 6

					echo "#### blastN on nt db - success ####"

					blastx -query $FILEPATH \
							-db $4 \
							-out ${FILENAME}_nt_blastx.txt\
							-num_threads 8 \
							-evalue 0.01 \
							-outfmt 6
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