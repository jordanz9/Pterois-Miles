# -*- coding: utf-8 -*-
'''
@author: jordanz

This python scripts saves the target sequence "contig_192_pilon" from the lionfish genome. Then extracts the dna sequences
based on the coordinates of "toxin_genes.txt" and writes the extracted sequences to "toxin_sequences.fasta".

'''


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


 
# Is used to check and verify the sequence IDs present in the FASTA file
fasta_file = "lionfish_assembly.fasta"
with open(fasta_file, "r") as file:
    for record in SeqIO.parse(file, "fasta"):
        seq_id = record.id
        sequence = str(record.seq)
        print(seq_id)


# Extracts and saves the target sequence ID contig_192_pilon
fasta_file = "lionfish_assembly.fasta"
target_seq_id = "contig_192_pilon"
output_file = "contig_192_pilon.fasta"


# Find the target sequence and save it to the output file
with open(fasta_file, "r") as file:
    genome_sequences = {}
    for record in SeqIO.parse(file, "fasta"):
        seq_id = record.id
        sequence = str(record.seq)
        genome_sequences[seq_id] = sequence

if target_seq_id in genome_sequences:
    sequence = genome_sequences[target_seq_id]
    with open(output_file, "w") as outfile:
        outfile.write(f">{target_seq_id}\n{sequence}\n")
    print(f"Sequence '{target_seq_id}' saved to {output_file}")
else:
    print(f"Target sequence {target_seq_id} not found in the genome.")

#%%
    
# Extracts dna sequences based on "best" hits' coordinates
contig_file = "contig_192_pilon.fasta"
toxins_file = "toxin_genes.txt"
output_file = "toxin_sequences.fasta"

sequences = SeqIO.to_dict(SeqIO.parse(contig_file, "fasta"))


toxins_coordinates = []
with open(toxins_file, "r") as file:
    header = next(file)  
    for line in file:
        line = line.strip()
        if line:
            columns = line.split("\t")
            if len(columns) >= 3: 
                try:
                    seq_id = columns[0]
                    start = int(columns[1])
                    end = int(columns[2])
                    if start > end: 
                        start, end = end, start
                    toxins_coordinates.append((seq_id, start, end))
                except ValueError:
                    continue  


with open(output_file, "w") as file:
    for seq_id, start, end in toxins_coordinates:
        if seq_id in sequences:
            sequence = sequences[seq_id].seq[start - 1:end]
            if start > end:
                sequence = sequence.reverse_complement()
                header = f" : {end}_{start}"
            else:
                header = f" : {start}_{end}"
            
            extracted_seq_record = SeqRecord(sequence, id=seq_id, description=header)
            SeqIO.write(extracted_seq_record, file, "fasta")

print("Extraction complete. Sequences written to:", output_file)
