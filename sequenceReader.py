import sys
from Bio import SeqIO

"""
	File Name: sequenceReader.py
	Author: Sidney Oderberg
	Description: This file is for reading DNA sequences
"""

# Read all of the sequences from a given .fasta file and returns
# them in a list
def read(file):

	sequences = []
	for seq_record in SeqIO.parse(file, "fasta"):
		sequences.append(seq_record.seq)

	return sequences




