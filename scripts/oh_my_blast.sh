#!/usr/bin/bash -i

# $2 contigs to blast

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
					# $4 : vrl_genbank.info

					blastn -query $2.fasta \
						  	-db $3 \
						  	-out $2.blastn.txt\
						  	-num_threads 8 \
						  	-evalue 0.01 \
						  	-outfmt 6

					echo "#### blastN - success ####"

					blastx -query $2.fasta \
							-db $3_prot \
							-out $2.blastx.txt\
							-num_threads 8 \
							-evalue 0.01 \
							-outfmt 6
					echo "#### blastX - success ####"

					echo "#### Virus identification ####"

					if [[ -s $2.blastn.txt ]]
					then 
						python3 ./blastn_virus_identify.py $2.blastn.txt $4 $2.blastn.taxon
					fi

					if [[ -s $2.blastx.txt ]]
					then 
						python3 ./blastx_virus_identify.py $2.blastx.txt $4 $2.blastx.taxon
					fi

					echo "#### DONE ####"
					;; 
				n|N) echo "Cancel ";; 
				*) echo "Cancel" ;; 
			esac
			shift
		;;
		
	-nrnt)
		echo "You're going to blast{n,x} your contigs over the nr/nt database" ;;
	esac
	shift
done

conda deactivate
echo "paleogenomic conda env OFF"