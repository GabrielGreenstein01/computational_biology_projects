import utils
import numpy as np
import pandas as pd
import os

class FragmentSet(object):
	def __init__(self, fragfile, rmsdfile):
		"""
		This class contains the fragment library for the input protein. It must do the following:
		- Read in fragment file and parse fragments at each position. Fragment files are of the form <protein>_<frag_len>mers.frag
		- Read in RMSD file containing pre-calculated RMSD to native structure for each fragment at each position.
		- Based on fragments and their corresponding RMSDs, rank fragments at each position by RMSD
		
		"""
		# Will contain list of (phi,psi) tuples for each frament in every position;
		# Each column will be a position and the i-th row will be a list of (phi,psi) tuples for 
		# the i-th fragment in that position
		self.fragDataSet = pd.DataFrame()

		# Will contain a sublist for each position, containing the indices
		# of fragments with the lowest to highest RMSD
		self.rmsdSortedFragments = []

		# Read lines from fragment file
		lines = open(fragfile, 'r').readlines()

		# Get position ranges by finding indices containing "position"
		indices = [i for i, s in enumerate(lines) if 'position' in s]
		indices.append(len(lines))

		positions = []

		# Group the fragments by position
		for i in range(0,len(indices)-1):
			positions.append(lines[indices[i]+1:indices[i+1]])

		# Keep track of current position; use to name column in self.fragDataSet
		curr_position = 0

		# Iterate through each position
		for position in positions:

			# Add to our tally
			curr_position += 1

			# Define "position_n_mers" which will contain the list of n-mers 
			# for each fragment for that position
			position_n_mers = []

			# n_mer of current fragment
			curr_n_mer = []

			# Iterate through all fragments in the position
			for i in range(1,len(position)):

				# Clean up string
				line = position[i].strip()

				if len(line) != 0:

					# Split line into array by spaces
					line = line.split()

					# Obtain phi, psi values
					phi = float(line[5])
					psi = float(line[6])

					# Add (phi,psi) tuple to current n-mer
					curr_n_mer.append(tuple((phi,psi)))

				# If len(line) == 0, we move onto the next fragment
				if len(line) == 0:

					# Add n-mer to position_column
					position_n_mers.append(curr_n_mer)

					# Re-initialize n-mer to start building the next n-mer
					curr_n_mer = []

			# Add position_column to dataframe
			df = pd.DataFrame({str(curr_position): position_n_mers})
			self.fragDataSet = pd.concat([self.fragDataSet, df], axis=1) 

			# Reninitialize to keep track of fragments for the next position
			position_n_mers = []


		# Read lines from RMSD file
		rmsd_lines = open(rmsdfile, 'r').readlines()

		# List will be composed of sublists containing the RMSD value for each fragment in a given position
		rmsd_positions = []

		# Determine number of positions 
		num_positions = int(rmsd_lines[-1::][0].strip().split()[0])

		# Initialize rmsd_positions to have the same number of sublists as positions
		for i in range(num_positions):
			rmsd_positions.append(list())

		# Iterate through each line of rmsd_file
		for i in range(len(rmsd_lines)):

			# Format line
			rmsd_line = rmsd_lines[i].strip()
			rmsd_line = rmsd_line.split()

			# Add rmsd value to appropriate position sublist
			rmsd_positions[int(rmsd_line[0]) - 1].append(rmsd_line[2])

		# Sort each position sub-list by arguments
		for i in range(num_positions):
			self.rmsdSortedFragments.append(np.argsort(rmsd_positions[i]))

	def get_lowRMS_fragments(self, pos, N):
		"""
		Returns the top-ranked fragments by RMSD at a defined position in the chain
		--------
		Params
			- pos (int): fragment position in chain (1-indexed)
			- N (int): number of fragments to return
		Returns
			- lowRMS_fragments (list): top N fragments at pos by RMSD. This should be a list of lists of (phi, psi) tuples. 
			  For example, a 3-mer fragment could be represented as the following: [(-60.892, 142.456), (-72.281, 128.933), (-132.337, -175.477)]
		"""

		# Get fragments at position from dataframe
		column = self.fragDataSet[str(pos)]

		# Sort column by top-N rmsd values
		argSort = self.rmsdSortedFragments[pos - 1][:N]

		# Return the (phi,psi) tuple list for each of these N fragments at the given position
		return column.loc[argSort].tolist()





