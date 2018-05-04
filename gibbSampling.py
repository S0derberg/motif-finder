import sys
import math
import time
import numpy as np
import utilities as utils
import sequenceReader as reader
import greedySearch as greedy

"""
	File Name: gibbSampling.py
	Author: Sidney Oderberg
	Description: This file implements the Gibb Sampling algorithm for motif finding
"""

LENGTH = 8
T = 500
INIT = None
BASE_DICT = {0: "A", 1: "C", 2: "G", 3: "T"}


# Find the likelihood of a given sequence based on a position weight matrix
def find_probability(seq, PWM):
	prob = 1
	for c in range(len(seq)):
		if seq[c] == "A":
			prob *= PWM[0,c]
		elif seq[c] == "C":
			prob *= PWM[1,c]
		elif seq[c] == "G":
			prob *= PWM[2,c]
		else:
			prob *= PWM[3,c]
	return prob


# Sample a new starting index using a sequence of probabilities representing the
# likelihood of each starting position
def sample(probs):
	rand = np.random.random()
	running_sum = 0
	for p in range(len(probs)):
		if rand > running_sum and rand < running_sum + probs[p]:
			return p
		else:
			running_sum += probs[p]
	return p


def gibbs_sampling(file, L, T, init=None, verbose=False):
	print("Loading DNA sequences.")
	sequenceList = reader.read(file)

	print("Beginning Gibbs Sampling Algorithm for Motif Finding")
	print("    Motif Length: {}".format(L))
	print("    Number of Samples: {}".format(T))
	print("    Number of Sequences: {}".format(len(sequenceList)))

	start_time = time.time()

	seq_1 = sequenceList[0]
	seq_2 = sequenceList[1]
	highest_IC = 0
	motif_list = []
	best_motif_list = []

	base_probs = utils.base_probabilities(sequenceList)

	if init == "greedy":
		motif_list, motif, IC, total_time = greedy.greedy_search(file, L)
	else:
		starting_positions = []
		for s in range(len(sequenceList)):
			index = np.random.random_integers(0, len(sequenceList[s])-L)
			starting_positions.append(index)
			motif_list.append(sequenceList[s][index:index+L])
	

	for t in range(T):
		sequence_z_index = np.random.random_integers(0,len(motif_list)-1)

		if sequence_z_index == 0:
			remaining_motifs = motif_list[1:]
		elif sequence_z_index == len(motif_list)-1:
			remaining_motifs = motif_list[:sequence_z_index]
		else:
			remaining_motifs = motif_list[:sequence_z_index] + motif_list[sequence_z_index+1:]

		PM = utils.profile_matrix(remaining_motifs, L)
		PWM = utils.position_weight_matrix(PM, len(remaining_motifs))

		sequence_z = sequenceList[sequence_z_index]
		probabilities = np.zeros(len(sequence_z)-L+1)
		for x in range(len(sequence_z)-L+1):
			Q_x = find_probability(sequence_z[x:x+L], PWM)
			probabilities[x] = Q_x

		probabilities = probabilities / np.sum(probabilities)

		new_index = sample(probabilities)
		new_motif = sequence_z[new_index:new_index+L]
		motif_list[sequence_z_index] = new_motif

		PM = utils.profile_matrix(motif_list, L)
		PWM = utils.position_weight_matrix(PM, len(motif_list))
		IC = utils.information_content(PWM, base_probs, L)

		if IC > highest_IC:
			highest_IC = IC
			best_motif_list = motif_list.copy()


	final_PM = utils.profile_matrix(best_motif_list, L)
	final_PWM = utils.position_weight_matrix(final_PM, len(sequenceList))
	final_IC = utils.information_content(final_PWM, base_probs, L)

	consensus_motif = utils.find_motif(PM)

	elapsed_time = time.time() - start_time

	print("\nGibbs Sampling Complete.")
	print("    Motif Found: " + consensus_motif)
	print("    Information Content: {}".format(final_IC))
	print("    Elapsed Time: {} seconds".format(elapsed_time))
	
	if verbose:
		for i in range(len(best_motif_list)):
			print(best_motif_list[i])
		print(final_PM)
		print(final_PWM)

	return consensus_motif, final_IC, elapsed_time


# Run the gibbs sampling algorithm h on the DNA sequence list
def main(args):

	motif, IC, time = gibbs_sampling("data/hm01r.fasta", LENGTH, T, init=INIT)
	

if __name__=="__main__":
	main(sys.argv)