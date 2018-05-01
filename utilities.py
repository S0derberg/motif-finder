import sys
import math
import numpy as np

"""
	File Name: utilities.py
	Author: Sidney Oderberg
	Description: This file consists of helper functions used by the main
				 algorithm modules.
"""


BASE_DICT = {0: "A", 1: "C", 2: "G", 3: "T"}

# Find the base probabilities of the nucleotides by looking at their frequences
# over all the sequences currently being looked at.
def base_probabilities(sequenceList):
	a_count = 0
	c_count = 0
	g_count = 0
	t_count = 0
	total = 0

	for seq in sequenceList:
		for s in range(len(seq)):
			total += 1
			if seq[s] == "A":
				a_count += 1
			elif seq[s] == "C":
				c_count += 1
			elif seq[s] == "G":
				g_count += 1
			else:
				t_count += 1

	return [a_count/total, c_count/total, g_count/total, t_count/total]

# From a list of sequences, make a profile matrix
def profile_matrix(sequences, L):
	PM = np.zeros((4,L))

	for i in range(L):
		a_count = 0
		c_count = 0
		g_count = 0
		t_count = 0

		for j in range(len(sequences)):
			if sequences[j][i] == "A":
				a_count += 1
			elif sequences[j][i] == "C":
				c_count += 1
			elif sequences[j][i] == "G":
				g_count += 1
			else:
				t_count += 1

		PM[0,i] = a_count
		PM[1,i] = c_count
		PM[2,i] = g_count
		PM[3,i] = t_count

	return PM


# From a Profile Matrix, make a Position Weight Matrix
def position_weight_matrix(PM, total):
	return PM / total


# From a Position Weight Matrix, find the Information Content. Higher is better
def information_content(PWM, base_probs, L):
	info_content = 0
	for k in range(L):
		for b in range(len(base_probs)):

			w_bk = PWM[b,k]
			q_b = base_probs[b]

			if w_bk == 0:
				continue

			info_content += w_bk * math.log((w_bk / q_b), 2)

	return info_content


# From a Profile Matrix, find the Consensus motif, a sequence where each
# nucleotide is the most likely at each position
def find_motif(PM):
	motif = ""
	for s in range(PM.shape[1]):
		highest = 0
		for b in range(4):
			if PM[b,s] > highest:
				highest = PM[b,s]
				base = BASE_DICT[b]

		motif += base

	return motif