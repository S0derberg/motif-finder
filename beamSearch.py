import sys
import math
import time
import numpy as np
import utilities as utils
import sequenceReader as reader

"""
	File Name: beamSearch.py
	Author: Sidney Oderberg
	Description: This file implements the Beam Search algorithm for motif finding
"""

LENGTH = 8
K = 10
BASE_DICT = {0: "A", 1: "C", 2: "G", 3: "T"}

# A data structure to keep track of the best K searches
class K_Best:
	def __init__(self, k):
		self.k = k
		self.data = {}
		self.worst = -1

		for i in range(k):
			self.data[i] = (-1, None)


	def add(self, value, data):
		if value <= self.worst:
			return

		node = (value, data)
		self.data[self.k-1] = node

		i = self.k-2
		while i >= 0:
			if value > self.data[i][0]:
				self.data[i+1] = self.data[i]
				self.data[i] = node
				i -= 1
			else:
				break
		self.worst = self.data[self.k-1][0]

	def get_data(self):
		return self.data

# Run the beam search on the DNA sequence list
def beam_search(file, L, K, verbose=False):
	print("Loading DNA sequences.")
	sequenceList = reader.read(file)

	print("Beginning Beam Search Algorithm for Motif Finding")
	print("    Motif Length: {}".format(L))
	print("    Number of Sequences: {}".format(len(sequenceList)))
	print("    K: {}".format(K))

	start_time = time.time()

	seq_1 = sequenceList[0]
	seq_2 = sequenceList[1]
	highest_IC = 0
	best_motif_list = []
	best_k = K_Best(K)

	base_probs = utils.base_probabilities(sequenceList)

	for i in range(len(seq_1)-L+1):
		motif_1 = seq_1[i:i+L]

		for j in range(len(seq_2)-L+1):
			motif_2 = seq_2[j:j+L]

			PM = utils.profile_matrix([motif_1, motif_2], L)
			PWM = utils.position_weight_matrix(PM, 2)
			IC = utils.information_content(PWM, base_probs, L)
			
			best_k.add(IC, [motif_1, motif_2])


	prev_best_k = best_k
	
	for k in range(len(sequenceList)-2):
		curr_seq = sequenceList[k+2]
		curr_best_k = K_Best(K)

		for p in range(len(curr_seq)-L+1):
			motif_k = curr_seq[p:p+L]

			for t in range(K):
				motif_list = prev_best_k.get_data()[t][1] + [motif_k]

				PM = utils.profile_matrix(motif_list, L)
				PWM = utils.position_weight_matrix(PM, k+3)
				IC = utils.information_content(PWM, base_probs, L)

				curr_best_k.add(IC, motif_list)

		prev_best_k = curr_best_k


	best_motif_list = prev_best_k.get_data()[0][1]

	final_PM = utils.profile_matrix(best_motif_list, L)
	final_PWM = utils.position_weight_matrix(PM, len(sequenceList))
	final_IC = utils.information_content(PWM, base_probs, L)

	consensus_motif = utils.find_motif(PM)

	elapsed_time = time.time() - start_time

	print("\nBeam Search Complete.")
	print("    Motif Found: " + consensus_motif)
	print("    Information Content: {}".format(final_IC))
	print("    Elapsed Time: {} seconds".format(elapsed_time))
	
	if verbose:
		for i in range(len(best_motif_list)):
			print(best_motif_list[i])
		print(final_PM)
		print(final_PWM)

	return consensus_motif, final_IC, elapsed_time


# Call the helper function to run the beam search
def main(args):

	motif, IC, time = beam_search("data/hm01r.fasta", LENGTH, K)
	


if __name__=="__main__":
	main(sys.argv)


