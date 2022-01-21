#!/bin/bash -i


mkdir seqref
cd seqref

conda activate paleogenomic

# get annotations
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.gpff.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.protein.gpff.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.protein.gpff.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.4.protein.gpff.gz
gunzip viral.1.protein.gpff.gz
gunzip viral.2.protein.gpff.gz
gunzip viral.3.protein.gpff.gz
gunzip viral.4.protein.gpff.gz

# download and create BLAST databases for viral refseq protein
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz
gunzip viral.1.protein.faa.gz
makeblastdb -dbtype prot \
        -parse_seqids \
        -in viral.1.protein.faa \
        -out viral.1.protein \
        -title 'Virus Refseq Protein 1'
        
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.protein.faa.gz
gunzip viral.2.protein.faa.gz
makeblastdb -dbtype prot \
        -parse_seqids \
        -in viral.2.protein.faa \
        -out viral.2.protein \
        -title 'Virus Refseq Protein 2'
        
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.protein.faa.gz
gunzip viral.3.protein.faa.gz
makeblastdb -dbtype prot \
        -parse_seqids \
        -in viral.3.protein.faa \
        -out viral.3.protein \
        -title 'Virus Refseq Protein 3'

wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.4.protein.faa.gz
gunzip viral.4.protein.faa.gz
makeblastdb -dbtype prot \
        -parse_seqids \
        -in viral.4.protein.faa \
        -out viral.4.protein \
        -title 'Virus Refseq Protein 4'
        
        
        
# Create an alias for the 4 DB :
blastdb_aliastool -dblist "viral.1.protein viral.2.protein viral.3.protein viral.4.protein" \
		-dbtype prot \
		-out viral.protein.all \
		-title "Viral Protein All"

conda deactivate