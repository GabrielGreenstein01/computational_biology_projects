import sys

import pandas as pd
import numpy as np

import csv
import math

import matplotlib.pyplot as plt
from collections import OrderedDict

class GSEA(object):
	"""
	GSEA class runs a geneset enrichment analysis on a given dataset
	"""
	def __init__(self):

		self.expression_file = ""
		self.sample_file = ""
		self.kegg_file = ""

		self.expressions = pd.DataFrame()
		self.samples = pd.DataFrame()

		self.labels = []

		self.genes = []
		self.gene_sets_dict = dict()

	def load_data(self, expfile, sampfile, genesets):
		"""
		Take the file paths to the expression, sample and gene set data, read them in and store within the GSEA instance
		"""

		self.expression_file = expfile
		self.sample_file = sampfile
		self.kegg_file = genesets

		# Read in X and Y
		self.expressions = pd.read_csv(expfile, sep='\t')
		self.samples = pd.read_csv(sampfile, sep='\t', header=None)

		# Get list of genes in expression file
		self.genes = self.expressions['SYMBOL'].tolist()
		self.genes = np.array(self.genes)

		# Remove genes from expression
		self.expressions = self.expressions.drop(columns=['SYMBOL'])

		indices = []

		# Reorder samples in terms of expression columns (i.e. the i-th column corresponds to the i-th row)
		for col in self.expressions:
			index = self.samples.loc[self.samples.iloc[:,0] == col].index.values[0]
			indices.append(index)
		
		self.samples = self.samples.reindex(indices)

		# Read in KEGG file
		with open(self.kegg_file, 'r') as f:
			reader = csv.reader(f, delimiter='\t')
			for row in reader:
				self.gene_sets_dict[row[0]] = row[2:len(row)]

	def reorder_labels(self,randomize):
		"""
		Provides a reordering of the labels (if randomize = 1) or resets the labels to their default (i.e., correct) values (if randomize = 0)
		"""

		# Set labels to be default values
		self.labels = np.array(list(self.samples[1]))

		# Randomly permute correct labels 
		if randomize == 1:
			r = np.random.default_rng()
			self.labels = r.permutation(self.labels)

	def get_gene_rank_order(self, randomize = None):
		"""
		Return a list of all genes (as strings) ranked by their logFC between patient and control,
		with the gene with the highest logFC ordered first
		"""
		# Reorders or reverts (back to correct values) labels depending on input
		if randomize is not None:
			self.reorder_labels(randomize)
		else:
			self.reorder_labels(0)

		# Get healthy/sick columns from the labels
		healthy_indices = np.where(self.labels == 0)[0]
		sick_indices = np.where(self.labels == 1)[0]

		# print(len(healthy_indices))
		# print(len(sick_indices))
		# print(np.where(self.genes == "BMP4"))
		# print(self.genes[967])
		# print(self.expressions.iloc[967,:])
		# print(self.expressions.iloc[967 , healthy_indices])
		# print(self.expressions.iloc[967 , sick_indices])

		# Calculate the mean value for each gene across the healthy and sick columns
		healthy = (self.expressions.iloc[: , healthy_indices].sum(axis=1)) / len(healthy_indices)
		sick = (self.expressions.iloc[: , sick_indices].sum(axis=1)) / len(sick_indices)

		# Compute their difference
		change = sick - healthy 

		# print(change[967])
		# print(change.sort_values(ascending=False))

		# Sort genes with highest differential expression first
		indices = (-change).argsort()
		indices = list(indices)

		# Find what genes those correspond to
		ordered_genes = self.genes[indices]

		# change = change[indices]
		# print(ordered_genes[566])
		# print(change[566])

		return ordered_genes

	def get_enrichment_score(self, geneset,randomize = None):
		"""
		Return the enrichment score, a float correct to two decimal places for a given gene set, such as ‘KEGG_CITRATE_CYCLE_TCA_CYCLE’
		(which is the string for the gene set name corresponding to the gene set). This method should run get_gene_rank_order at some point
		to initiate enrichment calculation
		"""
		# Get gene rank order for random labels or for original labels (depends on input)
		if randomize is not None:
			ordered_genes = self.get_gene_rank_order(randomize)
		else:
			ordered_genes = self.get_gene_rank_order()

		# Find genes in given geneset
		gene_set = self.gene_sets_dict[geneset]

		N = len(self.genes)
		G = 0

		# Only consider genes in the expression AND KEGG file
		for i in self.genes:
			if i in gene_set:
				G += 1

		# Define up and down movements in random walk
		up = math.sqrt((N-G)/G)
		down = math.sqrt(G/(N-G))

		running_sum = []
		current_sum = 0

		# Compute brownian bridge
		for gene in ordered_genes:

			if gene in gene_set:
				current_sum += up
			else:
				current_sum -= down

			running_sum.append(current_sum)

		# Calculate ES score
		ES = round(max(running_sum),2)

		return ES

	def get_sig_sets(self, p):
		"""
		Return the list of significant gene sets (as strings), at a corrected threshold of p, by name. If no gene sets are significant,
		return an empty list. This method should run get_gene_rank_order and/or get_enrichment_score at some point to initiate enrichment
		calculations and then identify significant gene sets.
		"""
		trials = 100

		sig_sets = []

		# Iterate through every geneset
		for key in self.gene_sets_dict:

			# Get enrichment score for original labels
			original_ES = self.get_enrichment_score(key)

			trial_ES = []

			# Run 100 trials on randomly permuted labels
			for i in range(0,trials):

				# Get enrichment score for randomly permuted labels
				ES = self.get_enrichment_score(key,1)
				trial_ES.append(ES)

			# Calculate the empirical p-value
			probability = len([i for i in trial_ES if i >= original_ES]) / trials

			# Implement Bonferroni correction
			corrected_p = probability / len(self.gene_sets_dict)

			# Determine if set is significant
			if corrected_p <= p:
				sig_sets.append(key)

		return sig_sets

	def calc_geneset_ES(self):

		# max_ES = [0]
		# max_ES_gs = ['']
		genes = dict()

		for key in self.gene_sets_dict:

			genes_in_set = self.gene_sets_dict[key]

			for i in range(0,len(genes_in_set)):
				genes[genes_in_set[i]] = 0


		for key in self.gene_sets_dict:

			genes_in_set = self.gene_sets_dict[key]

			for i in range(0,len(genes_in_set)):
				# if genes_in_set[i] in self.genes:
				genes[genes_in_set[i]] += 1
		# print(genes)


		# print(genes)
		max_gene = [""]
		occurences = 0
		for key in genes:
			if genes[key] > occurences:
				max_gene.clear()
				occurences = genes[key]
				max_gene.append(key)
			elif genes[key] == occurences:
				max_gene.append(key)
			else:
				continue

		print(max_gene)

		# len_genesets = dict()
		# for key in self.gene_sets_dict:
		# 	length = len(self.gene_sets_dict[key])
		# 	len_genesets[key] = length


		alphabetized_dictionary = OrderedDict(sorted(genes.items()))
		# print(alphabetized_dictionary)
		plt.bar(list(alphabetized_dictionary.keys()), alphabetized_dictionary.values(), color='royalblue')
		plt.title("Frequency of Genes in Gene Sets")
		plt.xlabel("Gene in KEGG (ordered a-z)")
		plt.ylabel("Number of Gene Sets")
		plt.xticks([])
		# plt.tick_params(axis='y', which='minor', bottom=False)
		plt.show()

			# ES = self.get_enrichment_score(key)

			# if ES > max_ES[0]:
			# 	max_ES.clear()
			# 	max_ES_gs.clear()

			# 	max_ES.append(ES)
			# 	max_ES_gs.append(key)
			# elif ES == max_ES[0]:
			# 	max_ES.append(ES)
			# 	max_ES_gs.append(key)
			# else:
			# 	continue
		return 


def main():
    """
    Function is called when the file is called in the terminal
    """
    # check that the file is being properly used
    if (len(sys.argv) !=4):
        print("Please specify an input file and an output file as args.")
        return
      
    # input variables
    expression_file = sys.argv[1]
    sample_file = sys.argv[2]
    kegg_file = sys.argv[3]


    gsea = GSEA()
    gsea.load_data(expression_file, sample_file, kegg_file)
    # gsea.reorder_labels(0)
    # print(np.where(gsea.get_gene_rank_order() == "BMP4")[0])
    print(gsea.calc_geneset_ES())
    # print(gsea.get_enrichment_score("KEGG_CITRATE_CYCLE_TCA_CYCLE"))
    # gsea.get_sig_sets(0.5)

    print("")

if __name__=="__main__":
    main()