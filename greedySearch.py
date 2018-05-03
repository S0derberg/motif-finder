import sys
import math
import time
import numpy as np
import utilities as utils
import sequenceReader as reader

"""
	File Name: greedySearch.py
	Author: Sidney Oderberg
	Description: This file implements the Greedy Search algorithm for motif finding
"""

L = 8
BASE_DICT = {0: "A", 1: "C", 2: "G", 3: "T"}

# Run the greedy search on the DNA sequence list
def greedy_search(file, verbose=False):
	print("Loading DNA sequences.")
	sequenceList = reader.read(file)

	print("Beginning Greedy Algorithm for Motif Finding")
	print("    Motif Length: {}".format(L))
	print("    Number of Sequences: {}".format(len(sequenceList)))

	start_time = time.time()

	seq_1 = sequenceList[0]
	seq_2 = sequenceList[1]
	highest_IC = 0
	best_motif_list = []

	base_probs = utils.base_probabilities(sequenceList)

	for i in range(len(seq_1)-L+1):
		motif_1 = seq_1[i:i+L]

		for j in range(len(seq_2)-L+1):
			motif_2 = seq_2[j:j+L]

			PM = utils.profile_matrix([motif_1, motif_2], L)
			PWM = utils.position_weight_matrix(PM, 2)
			IC = utils.information_content(PWM, base_probs, L)
			
			if IC > highest_IC:
				highest_IC = IC
				best_motif_list = [motif_1, motif_2]

	
	for k in range(len(sequenceList)-2):
		curr_seq = sequenceList[k+2]
		highest_IC = 0
		best_motif = None

		for p in range(len(curr_seq)-L+1):
			motif_k = curr_seq[p:p+L]

			PM = utils.profile_matrix(best_motif_list + [motif_k], L)
			PWM = utils.position_weight_matrix(PM, k+3)
			IC = utils.information_content(PWM, base_probs, L)

			if IC > highest_IC:
				highest_IC = IC
				best_motif = motif_k

		best_motif_list.append(best_motif)


	final_PM = utils.profile_matrix(best_motif_list, L)
	final_PWM = utils.position_weight_matrix(PM, len(sequenceList))
	final_IC = utils.information_content(PWM, base_probs, L)

	consensus_motif = utils.find_motif(PM)

	elapsed_time = time.time() - start_time

	print("\nGreedy Search Complete.")
	print("    Motif Found: " + consensus_motif)
	print("    Information Content: {}".format(final_IC))
	print("    Elapsed Time: {} seconds".format(elapsed_time))
	
	if verbose:
		for i in range(len(best_motif_list)):
			print(best_motif_list[i])
		print(final_PM)
		print(final_PWM)

	return best_motif_list


# Call the function that runs the greedy search
def main(args):

	motif_list = greedy_search("data/hm01r.fasta")
	


if __name__=="__main__":
	main(sys.argv)