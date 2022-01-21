#!/bin/bash -i

mkdir seqref
cd seqref

conda activate paleogenomic

# get annotations
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.genomic.gbff.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.genomic.gbff.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.genomic.gbff.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.4.genomic.gbff.gz
gunzip viral.1.genomic.gbff.gz
gunzip viral.2.genomic.gbff.gz
gunzip viral.3.genomic.gbff.gz
gunzip viral.4.genomic.gbff.gz


# download and create BLAST databases for viral refseq genomic
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
gunzip viral.1.1.genomic.fna.gz
makeblastdb -dbtype nucl \
        -parse_seqids \
        -in viral.1.1.genomic.fna \
        -out viral.1.1.genomic \
        -title 'Virus Refseq Genomic 1'
        
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz
gunzip viral.1.1.genomic.fna.gz
makeblastdb -dbtype nucl \
        -parse_seqids \
        -in viral.2.1.genomic.fna \
        -out viral.2.1.genomic \
        -title 'Virus Refseq Genomic 2'
        
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.1.genomic.fna.gz
gunzip viral.1.1.genomic.fna.gz
makeblastdb -dbtype nucl \
        -parse_seqids \
        -in viral.3.1.genomic.fna \
        -out viral.3.1.genomic \
        -title 'Virus Refseq Genomic 3'
        
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.4.1.genomic.fna.gz
gunzip viral.1.1.genomic.fna.gz
makeblastdb -dbtype nucl \
        -parse_seqids \
        -in viral.4.1.genomic.fna \
        -out viral.4.1.genomic \
        -title 'Virus Refseq Genomic 4'

conda deactivate