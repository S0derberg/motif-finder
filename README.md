# CS 466 Final Project: Motif Finding

This is the GitHub repository for my CS 466: Bioinformatics Final Project at the University of Illinois, Urbana-Champaign.  For my project I chose to compare the motif-finding algorithms that we learned in class: Greedy Search, Beam Search, and Gibbs Sampling.

## How to Run the Code

First, clone the repository.  To run the code for an individual algorithm, navigate to the main (outermost) directory and in a terminal/command prompt run the command "python greedySearch.py".  That will run the greedy search algorithm with default motif length 8 and DNA file "hm01r.fasta".  Those values can be changed at the top of each method's individual file.  Similarly, run "python beamSearch.py" or "python gibbSampling.py" for the other methods.  Each of these will find a motif from the file of DNA sequences and report the motif, information content, and program execution time.

To run the experiments that I ran for the report, use the command "python experiments.py."  This will run 5 experiments that compare the results of the different methods when values of motif length, K for beam search, and T (number of samples) for Gibbs sampling vary.  This program will plot the results that I put in my report.

## More Information

An explanation of the motif-finding problem and the algorithms I used, as well as my results and discussion, can be found in my final report, which is also in this repository as a PDF document.  The data I used can be found in the data folder, and my result plots are in the results folder.
