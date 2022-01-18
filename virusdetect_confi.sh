#!/usr/bin/sh

wget https://github.com/kentnf/VirusDetect/archive/refs/tags/v1.7.tar.gz
tar -xzvf v1.7.tar.gz
rm -rf v1.7.tar.gz

cd VirusDetect-1.7/databases/
rm -rf ./*

wget bioinfo.bti.cornell.edu/ftp/program/VirusDetect/virus_database/v239/U100/plant_239_U100.tar.gz
tar -xzvf plant_239_U100.tar.gz --strip-components=1
rm -rf plant_239_U100.tar.gz

wget bioinfo.bti.cornell.edu/ftp/program/VirusDetect/virus_database/v239/vrl_genbank.info.gz
wget bioinfo.bti.cornell.edu/ftp/program/VirusDetect/virus_database/v239/vrl_idmapping.gz
wget https://raw.githubusercontent.com/neonicoo/paleovirology-5BiM-project/main/ehraharta_erecta_seq.fasta

# sudo dnf install perl cpanminus perl-GD expat-devel #if Fedora or CentOS distro
# cpanm XML::Parser
# cpanm Bio::Seq
# cpanm Bio::Perl
# cpanm Bio::DB::SeqFeature::Store
# cpanm Bio::Graphics
