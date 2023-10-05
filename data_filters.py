# -*- coding: utf-8 -*-
'''
@author: jordanz

This Python script processes tBLASTn results, filters data, and extracts non-overlapping gene coordinates using functions:
def merge_coordinates() and  def non_overlapping().

Inputs:
"blast_results.txt" file

Outputs:
"toxin_genes.txt" file: the final filtered and processed data, containing non-overlapping gene coordinates

'''

from _parser import parse_blast_results
import pandas as pd

# %%
# functions to be used

def merge_coordinates(group):
    coordinates = []
    for start, end in zip(group['qstart'], group['qend']):
        coordinates.append(f"{start}-{end}")
    return ', '.join(coordinates)


def non_overlapping(df, start, end):
    
    non_overlapping_entries = []
    prev_subend = -1
    
    for i, row in df.iterrows():
        # Check if the current substart is greater than or equal
        # to the previous subend
        if row[start] >= prev_subend:
            non_overlapping_entries.append(row)
            prev_subend = row[end]
            
    non_overlapping_entries = pd.DataFrame(non_overlapping_entries) 
    out = non_overlapping_entries.reset_index(drop=True)
    
    return out
# %%

# Parse "blast_results" and create a DataFrame
br = parse_blast_results("blast_results.txt")

# Keep rows with 'sseqID' equal to 'contig_192_pilon'
contig_df = br[br['sseqID'] == 'contig_192_pilon']
    

# %%
# Filtered hits criteria
min_bitscore = 60
min_per_ident = 50


filtered_contig_df = contig_df[
    (contig_df['per_ident'] >= min_per_ident) &
    (contig_df['bitscore'] >= min_bitscore)
]


# %%
# Select desired columns and save to another file

selected_columns = ['sseqID', 'per_ident', 'qstart',
                    'qend', 'substart', 'subend']
selected_df = filtered_contig_df[selected_columns]


# Swap 'substart' and 'subend' values if 'substart' > 'subend'
selected_df.loc[selected_df['substart'] > selected_df['subend'],
                ['substart', 'subend']] = selected_df.loc[selected_df['substart']
                                                          > selected_df['subend'],
                                                          ['subend', 'substart']].values


# Sort the DataFrame based on 'substart' column in ascending order
sorted_df = selected_df.sort_values(by='substart', ascending=True)

# %%

# Group by 'sseqID', 'substart', and 'subend' and apply the merge_coordinates function
merged_df = selected_df.groupby(['sseqID', 'substart', 'subend']).apply(merge_coordinates).reset_index(name='coordinates')


# %%
 
# Reset the index of the new DataFrame
non_overlapping_gene_data = non_overlapping(merged_df, 'substart', 'subend')

# Display the non-overlapping gene data
print(non_overlapping_gene_data)

non_overlapping_gene_data.to_csv('toxin_genes.txt', sep = '\t', index = False)
