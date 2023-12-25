import pandas as pd
import numpy as np
from itertools import combinations
from itertools import product


def calculate_tanimoto(drugs_to_vec,drugA, drugB):
	"""
	Function calculates the Tanimoto coefficients of two drugs
	Inputs: drug => fingerprint dictionary, Drug A, Drug B
	Returns: tanimoto coefficient
	"""

	setA = drugs_to_vec[drugA]
	setB = drugs_to_vec[drugB]

	a = len(setA)
	b = len(setB)
	c = len(setA.intersection(setB))

	Tc = c / (a + b - c)

	return Tc

def drugs_to_vec(drugs_file):
	"""
	Function generates drug => fingerprint dictionary using the drugs.csv file
	Inputs: drugs.csv filename
	Returns: drug => fingerprint dictionary
	"""

	# Read drugs_file into dataframe
	drugs = pd.read_csv(drugs_file)  

	# Format drugs dataframe
	drugs = drugs.drop('generic_name', axis=1)
	drugs['maccs'] = pd.Series(drugs['maccs']).str.split().apply(set)

	# Convert to dictionary for faster lookup
	drugs_to_vec = drugs.set_index('db_id')['maccs'].to_dict()

	# Return dictionary
	return drugs_to_vec

def targets_df(targets_file):
	"""
	Function returns targets dataframe from file
	"""

	return pd.read_csv(targets_file)

def generate_combos(setA, setB = None):
	"""
	Function generates all drug combinations;
	If one set is provided, it returns all combinations
	If two sets are provided, it returns all product combinations, i.e., (drugAi, drugBj)
	"""

	if setB is None:
		return list(combinations(setA, 2))
	else:
		return list(product(setA, setB))




