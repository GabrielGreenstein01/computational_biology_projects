import pandas as pd
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt

import chemoUtils

import sys

class Tanimoto(object):
	"""
	Object computes the Tanimoto coefficients for all drug combinations in drugscsv
	"""

	def __init__(self, drugs_to_vec, targets_df):
		"""
		To reduce redundant operations, we require drug => fingerprint dictionary and target dataframe as inputs
		We use this to create a drug => target(s) dictionary which we store within the object
		"""

		# Store drugs_to_vec and targets_df from input
		self.drug_to_vec = drugs_to_vec
		self.targets_df = targets_df

		# Initialize drug => target(s) dictionary
		self.drug_to_target = dict()

		# Get list of unique drugs appearing in drugs.csv
		self.unique_drugs = list(self.drug_to_vec.keys())

		# Fill in drug => target(s) dictionary
		for drug in self.unique_drugs:
			self.drug_to_target[str(drug)] = set(self.targets_df.loc[self.targets_df['db_id'] == str(drug)]['uniprot_accession'])

	def share_target(self, drugA, drugB):
		"""
		Function checks whether two drugs have any common proteins to which they bind
		Input: Drug A, Drug B
		Output: 1 (had common protein), 0 (no shared protein)
		"""
		# Get sets of proteins each drug binds to 
		setA = self.drug_to_target[str(drugA)]
		setB = self.drug_to_target[str(drugB)]

		# Check if the two drugs share any common proteins to which they bind
		if len(setA.intersection(setB)) == 0:
			return 0
		else:
			return 1  

	def compute(self):
		"""
		Function computes the Tanimoto coefficient for each combination of drugs;
		Returns dataframe ready to be exported to a csv file
		"""

		# Find all drug combinations to compute tanimoto coefficients
		drug_combos = chemoUtils.generate_combos(self.unique_drugs)

		# Initialize empty dictionary which we will use to generate dataframe
		dictionary_data = dict()

		# Iterate through all drug combos
		for i in range(len(drug_combos)):

			drug = drug_combos[i]

			drugA = drug[0]
			drugB = drug[1]

			# Compute T_c and determine whether or not the drugs share a target
			tanimoto_coefficient = chemoUtils.calculate_tanimoto(self.drug_to_vec,drugA,drugB)
			share_target = self.share_target(drugA,drugB)

			# Append to dictionary
			dictionary_data[i] = [str(drugA), str(drugB), round(tanimoto_coefficient,6), share_target]

		# Convert dictionary to dataframe
		export_data = pd.DataFrame.from_dict(dictionary_data,orient='index')

		# Return dataframe
		return export_data

	def export(self, export_data, export_file):
		"""
		Function exports dataframe to csv file
		Inputs: dataframe to be exported, export filename
		"""

		# Export data to csv
		export_data.to_csv(export_file,index=False, header=False)

	def generate_histograms(self, export_data):
		"""
		Function generates histograms with appropriate bin sizes and labels
		"""

		# Generate "all_tanimoto.png"
		all_tanimoto = export_data[2].tolist()
		plt.hist(x=all_tanimoto, bins=20)

		plt.title("gbg222 All")
		plt.xlabel("Tanimoto Coefficient")
		plt.ylabel("Frequency")
		# plt.show()
		plt.savefig('all_tanimoto.png')

		plt.clf()

		# Generate "shared_tanimoto.png"
		shared_tanimoto = export_data.loc[export_data[3] == 1][2].tolist()
		plt.hist(x=shared_tanimoto, bins=20)
		plt.title("gbg222 Shared")
		plt.xlabel("Tanimoto Coefficient")
		plt.ylabel("Frequency")
		# plt.show()
		plt.savefig('shared_tanimoto.png')

		plt.clf()

		# Generate "notshared_tanimoto.png"
		not_shared_tanimoto = export_data.loc[export_data[3] == 0][2].tolist()
		plt.hist(x=not_shared_tanimoto, bins=20)
		plt.title("gbg222 Not Shared")
		plt.xlabel("Tanimoto Coefficient")
		plt.ylabel("Frequency")
		# plt.show()
		plt.savefig('notshared_tanimoto.png')

def main():
	"""
	Function is called when the file is called in the terminal
	"""
	# check that the file is being properly used
	if (len(sys.argv) !=4):
		print("Please specify an input file and an output file as args.")
		return

	# input variables
	drugs_file = sys.argv[1]
	targets_file = sys.argv[2]
	output_file = sys.argv[3]

	# Get dictionary drug => fingerprint and target dataframe from utility file
	drugs_to_vec = chemoUtils.drugs_to_vec(drugs_file)
	targets_df = chemoUtils.targets_df(targets_file)

	# Run Tanimoto functions to compute and export the data; also generate histograms
	tanimoto = Tanimoto(drugs_to_vec, targets_df)
	data = tanimoto.compute()
	tanimoto.export(data, output_file)
	tanimoto.generate_histograms(data)

if __name__=="__main__":
    main()