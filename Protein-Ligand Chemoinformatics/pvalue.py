import pandas as pd
import numpy as np
from itertools import product
import argparse
import random

import chemoUtils

import sys

class ligandSetSimilarity(object):
	"""
	Object returns the bootstrap p-value of two given drugs
	"""

	def __init__(self, drugs_to_vec, targets_df, proteinA, proteinB, n, r):
		"""
		Initialize object with drug => fingerprint dictionary, targets dataframe, protein A,
		protein B, and the number of sampling iterations (n)
		"""

		# Store number of iterations and random seed generator
		self.n = n

		# Initialize random number generator
		random.seed(r)

		# Store targets dataframe and drugs => fingerprint dictionary 
		self.targets_df = targets_df
		self.drugs_to_vec = drugs_to_vec

		# Get set of ligands which bind to proteins A and B
		self.setA = set(self.targets_df.loc[self.targets_df['uniprot_accession'] == str(proteinA)]['db_id'])
		self.setB = set(self.targets_df.loc[self.targets_df['uniprot_accession'] == str(proteinB)]['db_id'])

		# Store length of set
		self.nA = len(self.setA)
		self.nB = len(self.setB)

	def compute_T_sum(self, setA, setB):
		"""
		Function computes T_summary of two sets
		Inputs: Set A, Set B
		"""

		# Get all combinations
		drug_combos = chemoUtils.generate_combos(setA, setB)

		# Store T_summary
		total = 0

		# For each combination, compute the Tanimoto coefficient;
		# only store (i.e., add to total) if value > 0.5
		for combo in drug_combos:

			drugA = combo[0]
			drugB = combo[1]

			T = chemoUtils.calculate_tanimoto(self.drugs_to_vec, drugA, drugB)

			if T > 0.5:
				total += T

		return total

	def sample(self):
		"""
		Function samples from drugs.csv; returns two random sets of size nA and nB
		"""

		# Initialize new ligand set
		newSetA = set()
		newSetB = set()

		# Get list all unique drugs; get length of this list
		# (use for upper bound of random int generator)
		drugs = list(self.drugs_to_vec.keys())
		num_drugs = len(drugs)

		# Fill newSetA with random drugs until it is the same size as nA
		while len(newSetA) < self.nA:
			intA = random.randint(0,num_drugs-1)	
			newSetA.add(drugs[intA])

		# Fill newSetB with random drugs until it is the same size as nB
		while len(newSetB) < self.nB:
			intB = random.randint(0,num_drugs-1)	
			newSetB.add(drugs[intB])		

		# Return new sets
		return newSetA, newSetB

	def compute_pval(self):
		"""
		Function computes the bootstrap p-value
		"""

		# Compute T_summary for the original ligand sets of protein A and B
		T = self.compute_T_sum(self.setA, self.setB)

		total = 0

		# Sample n times
		for i in range(self.n):

			# Get new sets with the same number of binding ligands
			newSetA, newSetB = self.sample()

			# Compute T_i, i.e. T_summary for these new sets
			Ti = self.compute_T_sum(newSetA, newSetB)

			# Add +1 to numerator of p_bootstrapp if T_i >= T_summary
			if Ti >= T:
				total += 1	

		# Compute p_boostrap
		p_bootstrap = total / self.n

		# Return p_bootstrap
		return p_bootstrap

def main():
	"""
	Function is called when the file is called in the terminal
	"""

	# Define Parser
	parser = argparse.ArgumentParser(prog = 'p-value calculator', description = 'Generate a bootstrap p-value for the comparison of two proteins')
	parser.add_argument('-n', default = 500, type = int, help = 'number of iterations')
	parser.add_argument('-r', default = 214, type = int, help = 'parameter that sets the state of the pseudo-random number generator in Python')
	parser.add_argument('drugs_file', type=str)
	parser.add_argument('targets_file', type=str)
	parser.add_argument('proteinA', type=str)
	parser.add_argument('proteinB', type=str)

	args = parser.parse_args()

	# Store variables from shell
	n = args.n
	r = args.r
	drugs_file = args.drugs_file
	targets_file = args.targets_file
	proteinA = args.proteinA
	proteinB = args.proteinB

	# Get drug => fingerprint dictionary and targets dataframe from chemoUtils file
	drugs_to_vec = chemoUtils.drugs_to_vec(drugs_file)
	targets_df = chemoUtils.targets_df(targets_file)

	# Run ligandSetSimilarity
	ligandSetSim = ligandSetSimilarity(drugs_to_vec, targets_df, proteinA, proteinB, n, r)
	print(ligandSetSim.compute_pval())

if __name__=="__main__":
	main()		

