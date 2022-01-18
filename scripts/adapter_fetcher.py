#!/usr/bin/env python
# coding: utf-8


from pathlib import Path
import pandas as pd
import os
import sys

raws = []
raw_dir = str(sys.argv[1])
for r in os.listdir(Path(raw_dir)):
    raws.append(r)

adapters_index = pd.read_table("./96_MID_BO1.txt", header = None)
adapters_index.columns = ["Index", "Adapter"]
adapters_index = adapters_index.drop([0])
adapters_index.index = range(len(adapters_index))


fastq_adapters = {}
for fastq in raws :
    fastq_adapters[fastq] = fastq_adapters.get(str(fastq), None)


for fastq in fastq_adapters.keys():
    index = fastq.split("_")[-1].split('.')[0]
    adapter = adapters_index.loc[adapters_index['Index'] == index]['Adapter'].values[0]
    fastq_adapters[fastq] = adapter
