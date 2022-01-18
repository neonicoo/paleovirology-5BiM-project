#!/home/nicolas/miniconda3/bin/python3

import pandas as pd
import sys 
import os

blast = pd.read_table(sys.argv[1], header = None)
blast.columns = ["query_seqid",
				"subject_seqid",
				"perc_ident",
				"length",
				"mismatch",
				"gapopen",
				"query_start",
				"query_end",
				"subject_start",
				"subject_end",
				"evalue",
				"bitscore"]


virus = pd.read_table(sys.argv[2], header = None)
virus.columns = ["seqid",
				"length",
				"family",
				"name",
				"realm",
				"order",
				"autority"]

output = sys.argv[3]

sseqid = list(blast["subject_seqid"].values)
sname = []
for s in sseqid:
    sname.append(virus["name"].loc[virus["seqid"] == s].values[0])

    
blast.insert(2,'subject_name', sname)
blast.to_csv("./"+str(output)+".txt", header=True, sep='\t')

print("######## DONE #######")