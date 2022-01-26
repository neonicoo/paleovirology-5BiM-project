#!/home/nicolas/miniconda3/bin/python3

import pandas as pd
import sys 
import os

blastx = pd.read_table(sys.argv[1], header = None)
blastx.columns = ["#query_seqid",
				"#subject_idmapping",
				"#perc_ident",
				"#length",
				"#mismatch",
				"#gapopen",
				"#query_start",
				"#query_end",
				"#subject_start",
				"#subject_end",
				"#evalue",
				"#bitscore"]


virus_mapp = pd.read_table(str(sys.argv[2])+"vrl_idmapping", header = None)
virus_mapp.columns = ["#seqid",
					"#idmapping"]

virus_gff = pd.read_table(str(sys.argv[2])+"vrl_genbank.info", header = None)
virus_gff.columns = ["#seqid",
					"#length",
					"#genus",
					"#name",
					"#realm",
					"#order",
					"#autority"]

output = sys.argv[3]

sidmapping = list(blastx["#subject_idmapping"].values)
sseqid = []
sname = []
sgenus = []
for s in sidmapping:
    seqid = virus_mapp["#seqid"].loc[virus_mapp["#idmapping"] == s].values[0]
    sseqid.append(seqid)
    sname.append(virus_gff["#name"].loc[virus_gff["#seqid"] == seqid].values[0])
    sgenus.append(virus_gff["#genus"].loc[virus_gff["#seqid"] == seqid].values[0])
    
blastx.insert(2,'#subject_species', sname)
blastx.insert(3,'#subject_genus', sgenus)
blastx.to_csv(str(output)+".txt", header=True, sep='\t')