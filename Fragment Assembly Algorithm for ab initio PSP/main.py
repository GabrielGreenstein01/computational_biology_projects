"""
This is the master file, which you should use to set up and run the simulations.
You may define functions or classes as necessary

For an input sequence, do the following (see project page for details):
	1. load sequence from fasta file and initialize protein into extended configuration
	2. Run a series of simulations (n=5 or 10):
		- Perform MCMC sampling with 9-mer fragments from kT=100 to kT=1 (assembly stage)
		- Perform MCMC sampling with 3-mer fragments from kT=1 to kT=0.1 (refinement stage)
		- Take best (lowest-energy) structure after refinement and perform energy minimization (see utils.relax)
		- Log energy and RMSD to native structure after minimization
	3. Visualize lowest-RMSD structure in PyMol

"""

from pyrosetta import *
from pyrosetta.rosetta import *
init(extra_options='-mute all -constant_seed')

from Protein import Protein
from FragmentSampler import MCMCSampler

import utils
import argparse
import numpy as np
import time
import os

def main():

	# Define Parser
	parser = argparse.ArgumentParser(prog = 'Predict Protein Structure', description = 'The following program predicts the structure of a protein')
	parser.add_argument('--fasta', default = None, help = '.fasta file containing sequence')
	parser.add_argument('--logdir', default = '.', help = 'directory to save all log files')
	parser.add_argument('--nsims', default = 1, help = 'number of simulations') 
	parser.add_argument('--nfrags', default = 3, help = 'number of fragments to sample from at each iteration')
	parser.add_argument('--anneal_rate', default = 0.999, help = 'temperature annealing parameter')

	args = parser.parse_args()

	# Create folder specified in shell
	folder_name = args.logdir

	new_directory = f"../{folder_name}"

	if not os.path.exists(new_directory):
		os.makedirs(new_directory)

	# Open fasta file (locally)
	fasta_file = open("../starter_data/" + args.fasta, 'r')

	# Open fasta file (Gradescope)
	# fasta_file = open(args.fasta, 'r')

	# Read file, get sequence, initialize protein
	lines = fasta_file.readlines()
	seq = lines[1].strip()

	protein = Protein(sequence = seq)

	# Get protein name, size of fragment sample space, and number of simulations
	protein_name = args.fasta.split(".")[0]	

	N = args.nfrags

	num_simulations = args.nsims

	# Store data after each simulation;
	# Will write to "simulation_summary.txt"
	data = []

	# Run through iterations
	for i in range(int(num_simulations)):

		# Make subfolder for each simulation
		simulation = i + 1

		sim_folder = f"../{folder_name}/sim_{simulation}"
		if not os.path.exists(sim_folder):
			os.makedirs(sim_folder)

		# Assemble protein structure from 9-mer sampling
		MCMC_9mer = MCMCSampler(protein_name,protein,9,N,sim_folder)
		protein_9mer = MCMC_9mer.simulate()

		# Refine output of 9-mer construction using 3-mers
		MCMC_3mer = MCMCSampler(protein_name,protein_9mer,3,N,sim_folder)
		protein_3mer = MCMC_3mer.simulate()

		# Perform energy minimization (locally)
		result = utils.relax(f'{sim_folder}/best.pdb', f'../starter_data/{protein_name}.pdb')

		# Perform energy minimization (Gradescope)
		# result = utils.relax(f'{sim_folder}/best.pdb', f'{protein_name}.pdb')

		# Append info to data
		data.append([str(simulation), str(result[2]), str(result[1])])


	# Write "simulation_summary.txt" using data (locally)
	filename = f'{new_directory}/simulation_summary.txt'

	# Write "simulation_summary.txt" using data (Gradescope)
	# filename = '/autograder/simulation_summary.txt'

	with open(filename, 'w') as f:
		for entry in data:
			f.write('\t'.join(entry[:]) + '\n')


if __name__=='__main__':
    main()

