import pandas as pd
import numpy as np
import argparse

import chemoUtils
from pvalue import ligandSetSimilarity

class Network(object):
	"""
	Object generates network_edgelist.text
	"""
	def __init__(self, drugs_file, targets_file, nodes_file):
		"""
		Initialize object with drugs_file, targets_file, and nodes_file
		"""

		# Get drug => fingerprint dictionary and targets dataframe
		self.drugs_to_vec = chemoUtils.drugs_to_vec(drugs_file)
		self.targets_df = chemoUtils.targets_df(targets_file)

		# Read in nodes_file and get list of protein names
		nodes_df = pd.read_csv(nodes_file, delimiter=',')
		self.proteins = nodes_df['uniprot_accession'].tolist()

	def generate_file(self,n,r):
		"""
		Function generates the data contained in network_edgelist.text 
		Inputs: # of sample iterations, random seed
		"""

		# Store data
		data = []

		# Generate all protein combos
		protein_combos = chemoUtils.generate_combos(self.proteins)

		# Calculate the bootstrap p-value for each combo;
		# store in data only if bootstrap p-val is significant (i.e. <= 0.05)
		for combo in protein_combos:
			proteinA = combo[0]
			proteinB = combo[1]
			ligandSetSim = ligandSetSimilarity(self.drugs_to_vec, self.targets_df, proteinA, proteinB, n, r)
			pval = ligandSetSim.compute_pval()

			if pval <= 0.05:
				data.append(str(proteinA) + ' ' + str(proteinB))

		# Return data
		return data

	def export(self, data):
		"""
		Function exports data to .txt file with desired format
		Input: data to export
		"""
		# Write data to file
		with open('network_edgelist.txt', 'w') as f:
		    for entry in data:
		        f.write(entry)
		        f.write('\n')

def main():
	"""
	Function is called when the file is called in the terminal
	"""
	# check that the file is being properly used

	# Define Parser
	parser = argparse.ArgumentParser(prog = 'p-value calculator', description = 'Generate a bootstrap p-value for the comparison of two proteins')
	parser.add_argument('-n', default = 500, type = int, help = 'number of iterations')
	parser.add_argument('-r', default = 214, type = int, help = 'parameter that sets the state of the pseudo-random number generator in Python')
	parser.add_argument('drugs_file', type=str)
	parser.add_argument('targets_file', type=str)
	parser.add_argument('protein_nodes', type=str)

	args = parser.parse_args()

	# Store variables from shell
	n = args.n
	r = args.r
	drugs_file = args.drugs_file
	targets_file = args.targets_file
	nodes_file = args.protein_nodes

	# Run network object, generate data, and export to .txt with correct formatting
	network = Network(drugs_file, targets_file, nodes_file)
	data = network.generate_file(n, r)
	network.export(data)

if __name__=="__main__":
	main()		

