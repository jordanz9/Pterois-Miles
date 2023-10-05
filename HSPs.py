# -*- coding: utf-8 -*-
'''
@author: jordanz

This Python script creates a scatter plot of High-Scoring Segment Pairs (HSPs) from BLAST results.
Separate those alignments with conting_192 from the others.

'''


import pandas as pd
import matplotlib.pyplot as plt

# Read the blast results into a DataFrame
df = pd.read_csv("blast_results.txt", sep=r"\s+")

# Define the criteria for filtering
criteria = (df['per_ident'] > 50) & (df['bitscore'] > 60) & (df['sseqID'] == 'contig_192_pilon')

# Create a scatter plot with color differentiation
plt.figure(figsize=(6, 4))

# Plot all matches in blue
plt.scatter(df[~criteria]['per_ident'], df[~criteria]['bitscore'], s=7, 
            marker='o', alpha=0.7, color='g', label='Rest HSPs')

# Plot specific matches in red
plt.scatter(df[criteria]['per_ident'], df[criteria]['bitscore'], s=7,
            marker='o', alpha=0.7, color='coral', label='192_pilon HSPs')

plt.xlabel('Percentage Identity (%)')
plt.ylabel('Bit Score')
plt.title('Scatter Plot of HSPs')
plt.legend()
plt.grid(True)
plt.tight_layout()


# Save the plot with high DPI (e.g., 300)
plt.savefig('HSP.png', dpi=300)


plt.show()
