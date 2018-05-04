import sys
import numpy as np
import matplotlib.pyplot as plt
import greedySearch as greedy
import beamSearch as beam
import gibbSampling as gibbs

"""
	File Name: experiments.py
	Author: Sidney Oderberg
	Description: This file is what I used to run experiments with my different search
				 implementations.
"""


# Run experiments using the different algorithms
def main(args):


	Ls = [6,8,14,26,33]
	Ks = [5,10,15,20,25]
	Ts = [10,50,100,500,1000]


	ICs = []
	times = []
	for length in Ls:
		motif_list, motif, IC, time = greedy.greedy_search("data/hm01r.fasta", length)
		ICs.append(IC)
		times.append(time)

	plt.scatter(Ls, ICs)
	plt.title("Greedy Search: Information Content vs. Motif Length")
	plt.xlabel("Motif Length (Number of Characters)")
	plt.ylabel("Information Content")
	plt.show()

	plt.scatter(Ls, times, c="g")
	plt.title("Greedy Search: Program Time vs. Motif Length")
	plt.xlabel("Motif Length (Number of Characters)")
	plt.ylabel("Execution Time (Seconds)")
	plt.show()

	ICs = []
	times = []
	for length in Ls:
		motif, IC, time = beam.beam_search("data/hm01r.fasta", length, Ks[1])
		ICs.append(IC)
		times.append(time)

	plt.scatter(Ls, ICs)
	plt.title("Beam Search: Information Content vs. Motif Length, K=10")
	plt.xlabel("Motif Length (Number of Characters)")
	plt.ylabel("Information Content")
	plt.show()

	plt.scatter(Ls, times, c="g")
	plt.title("Beam Search: Program Time vs. Motif Length, K=10")
	plt.xlabel("Motif Length (Number of Characters)")
	plt.ylabel("Execution Time (Seconds)")
	plt.show()

	ICs = []
	times = []
	for k in Ks:
		motif, IC, time = beam.beam_search("data/hm01r.fasta", Ls[1], k)
		ICs.append(IC)
		times.append(time)

	plt.scatter(Ks, ICs)
	plt.title("Beam Search: Information Content vs. K, L=8")
	plt.xlabel("K")
	plt.ylabel("Information Content")
	plt.show()

	plt.scatter(Ks, times, c="g")
	plt.title("Beam Search: Program Time vs. K, L=8")
	plt.xlabel("K")
	plt.ylabel("Execution Time (Seconds)")
	plt.show()

	ICs = []
	times = []
	for length in Ls:
		motif, IC, time = gibbs.gibbs_sampling("data/hm01r.fasta", length, Ts[3])
		ICs.append(IC)
		times.append(time)

	plt.scatter(Ls, ICs)
	plt.title("Gibbs Sampling: Information Content vs. Motif Length, T=500")
	plt.xlabel("Motif Length (Number of Characters)")
	plt.ylabel("Information Content")
	plt.show()

	plt.scatter(Ls, times, c="g")
	plt.title("Gibbs Sampling: Program Time vs. Motif Length, T=500")
	plt.xlabel("Motif Length (Number of Characters)")
	plt.ylabel("Execution Time (Seconds)")
	plt.show()

	ICs = []
	times = []
	slowIC = []
	slowTimes = []
	for t in Ts:
		motif, IC, time = gibbs.gibbs_sampling("data/hm01r.fasta", Ls[1], t)
		ICs.append(IC)
		times.append(time)

		motif2, IC2, time2 = gibbs.gibbs_sampling("data/hm01r.fasta", Ls[1], t, init="greedy")
		slowIC.append(IC2)
		slowTimes.append(time2)

	plt.scatter(Ts, ICs)
	plt.scatter(Ts, slowIC, c="r")
	plt.title("Gibbs Sampling: Information Content vs. T, L=8, Both Starting Methods")
	plt.xlabel("T (Number of Samples)")
	plt.ylabel("Information Content")
	plt.show()

	plt.scatter(Ts, times, c="g")
	plt.scatter(Ts, slowTimes, c="r")
	plt.title("Gibbs Sampling: Program Time vs. T, L=8, Both Starting Methods")
	plt.xlabel("T (Number of Samples)")
	plt.ylabel("Execution Time (Seconds)")
	plt.show()






if __name__=="__main__":
	main(sys.argv)

