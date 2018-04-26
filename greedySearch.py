import sys
import sequenceReader as reader

"""
	File Name: greedySearch.py
	Author: Sidney Oderberg
	Description: This file implements the Greedy Search algorithm for motif finding
"""

L = 10

# Load in a list of sequences from one of the files in the data folder.
def main(args):

	sequenceList = reader.read("data/hm01r.fasta")

	for seq in sequenceList:
		print(seq)
		print()





if __name__=="__main__":
	main(sys.argv)