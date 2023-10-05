# -*- coding: utf-8 -*-
'''
@author: jordanz

This Python script reads 'blast_results.txt' and returns the parsed data as a Pandas DataFrame.

'''

import pandas as pd
import re

# Function to parse "blast_results" and return a DataFrame
def parse_blast_results(blast_results_file):
    data = {
        'qseqID': [],
        'sseqID': [],
        'per_ident': [],
        'length': [],
        'mismatch': [],
        'gap': [],
        'qstart': [],
        'qend': [],
        'substart': [],
        'subend': [],
        'evalue': [],
        'bitscore': []
    }

    with open(blast_results_file, "r") as f:
        # Skip header line if present
        header = f.readline().strip()
        if "per_ident" not in header:
            f.seek(0)  # Reset file pointer if no header found

        for line in f:
            if line.startswith("#"):
                # Skip comment lines
                continue
            fields = re.split(r'\s+', line.strip())
            if len(fields) != 12:
                # Skip lines that don't have the expected number of columns
                continue

            data['qseqID'].append(fields[0])
            data['sseqID'].append(fields[1])
            data['per_ident'].append(float(fields[2]))
            data['length'].append(int(fields[3]))
            data['mismatch'].append(int(fields[4]))
            data['gap'].append(int(fields[5]))
            data['qstart'].append(int(fields[6]))
            data['qend'].append(int(fields[7]))
            data['substart'].append(int(fields[8]))
            data['subend'].append(int(fields[9]))
            data['evalue'].append(float(fields[10]))
            data['bitscore'].append(float(fields[11]))

    return pd.DataFrame(data)

blast_results_df = parse_blast_results("blast_results.txt")
