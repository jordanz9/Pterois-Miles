# -*- coding: utf-8 -*-
'''
@author: jordanz

This Python script reads the tBLASTn result file ('blast_results.txt') and produces three plots:
(i) Identity vs sseqID, (ii) Frequency of sseqID, (iii) Bitscore vs sseqID

For these plots, a threshold frequency was calculated to filter out less frequent sseqID values.
Threshold frequency is set at 4% of the highest frequency in the dataset.

'''

import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("blast_results.txt", sep=r"\s+")

# Calculate the frequency of each sseqID
sseqid_counts = df['sseqID'].value_counts()

# Calculate the threshold frequency (4% of the largest frequency)
threshold_frequency = 0.04 * sseqid_counts.max()

# Filter out sseqID values that are smaller than the threshold frequency
filtered_sseqid_counts = sseqid_counts[sseqid_counts >= threshold_frequency]

plt.figure(figsize=(16, 8))

# Identity
plt.figure(1)
filtered_df = df[df['sseqID'].isin(filtered_sseqid_counts.index)]
plt.scatter(filtered_df['per_ident'], filtered_df['sseqID'], marker='o',
            s = 7, color = 'g', alpha=0.7)
plt.xlabel('Identity (%)')
plt.ylabel('sseqID')
plt.title('Identity vs sseqID (Filtered)')
plt.grid(False)
plt.savefig('id.png', dpi=300)
plt.show()


# Frequency
plt.figure(2)
plt.barh(filtered_sseqid_counts.index, filtered_sseqid_counts.values,
         color='g', edgecolor='k')
plt.xlabel('Frequency')
plt.ylabel('sseqID')
plt.title('Frequency of sseqID (Filtered)')
plt.grid(False)
plt.savefig('freq.png', dpi=300)
plt.show()

# Bitscore
plt.figure(3)
filtered_df = df[df['sseqID'].isin(filtered_sseqid_counts.index)]
plt.scatter(filtered_df['bitscore'], filtered_df['sseqID'], marker='o',
            s = 7, color = 'g', alpha=0.7)
plt.xlabel('Bitscore')
plt.ylabel('sseqID')
plt.title('Bitscore vs sseqID (Filtered)')
plt.grid(False)

# Adjust the layout
plt.tight_layout()
plt.savefig('bitscore.png', dpi=300)
plt.show()
