"""
This file contains the main fragment sampling class, which performs a Monte Carlo simulated annealing procedure to fold a protein.
"""

from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.scoring import *
init(extra_options='-mute all  -constant_seed')

import numpy as np
import utils
from Protein import Protein
from FragmentSet import FragmentSet

import os
import random


class MCMCSampler(object):
    def __init__(self, protein_name, protein, frag_length, N, sim_folder):
        """
        Initialize necessary variables
        The score function is given to you (Rosetta centroid score function)
        **** NOTE: uncommented lines works locally but not on Gradescope
        --------
        Params:
            - protein_name (str): name of the protein we are predicting; use to name and import necessary files
            - protein (Protein): protein object we are sampling; allows us to refine 9-mer with 3-mer fragments
            - frag_length (int): length of fragments (9-mer or 3-mer); used to name and import necessary files
            - N (int): Size of fragment sample space (we only consider top N fragments)
            - sim_folder (str): string for where we save files
        """
        # Initialize scoring function (provided)
        self.scorefxn = create_score_function('score3')

        self.frag_length = frag_length
        self.sim_folder = sim_folder

        self.protein = protein

        # Get filenames (locally)
        frag_filename = "../starter_data/" + protein_name + f"_{frag_length}mers.frag"
        rmsd_filename = "../starter_data/" + protein_name + f"_{frag_length}mers.rmsd"

        # Get filenames (Gradescope)
        # frag_filename = protein_name + f"_{frag_length}mers.frag"
        # rmsd_filename = protein_name + f"_{frag_length}mers.rmsd"

        # Create FragmentSet
        self.fragments = FragmentSet(frag_filename, rmsd_filename)

        self.N = int(N)

        # Get column names from fragments dataset; use to find range of positions to sample from
        fragment_col_names = self.fragments.fragDataSet.columns.to_numpy()
        a = int(fragment_col_names[0])
        b = int(fragment_col_names[-1::])

        self.pos_range = [a,b]


        # Define dictionary to store top N fragments
        # Only gets top N fragments once to avoid redundant operations
        self.N_fragments = dict()

        # Get top N fragments at each position and store in dictionary
        for position in fragment_col_names.astype(int):

            self.N_fragments[position] = self.fragments.get_lowRMS_fragments(position, self.N)

        # Define annealing temperature parameters based on fragment length
        self.T_start = 0
        self.T_end = 0

        if frag_length == 9:
            self.T_start = 100
            self.T_end = 1
        if frag_length == 3:
            self.T_start = 1
            self.T_end = 0.1

        self.curr_T = self.T_start
        self.annealing_rate = 0.999

        # Define a set of sampled candidates at a given position at the same step
        self.sampled_candidates = set()

    def compute_energy(self, protein):
        """
        Compute energy of protein.
        --------
        Params:
            - protein (Protein object): protein to score
        Return:
            - energy of conformation (float)
        """        
        return utils.score_pose(protein, self.scorefxn)

    def perturb_fragment(self, position):
        """
        Sample from possible fragments for a position, and replace torsion angles of that fragment in the protein.
        ---------
        Params:
            - position (int): position in protein where we want to insert fragment
        Returns:
            - *copied* protein containing insertion
        """

        # Generate new seed and sample out of N fragments
        random.seed()
        sample_num = random.randint(0,self.N-1)

        # Keep track of candidate sampled; will reset once we move onto the next step
        # whether we accepted the insertion or rejected it and exhausted all other possible fragments
        self.sampled_candidates.add(sample_num)

        # Get fragment from sample space
        sample_fragment = self.N_fragments[position][sample_num]

        # Copy protein by allocating new memory
        copy_protein = Protein(pose = self.protein.pose)

        # Set (phi,psi) values of copied protein at positions to be equal to the fragment's (phi,psi) values
        for i in range(len(sample_fragment)):
            phi = sample_fragment[i][0]
            psi = sample_fragment[i][1]

            copy_protein.set_torsion(position + i, phi, psi)

        # Return copied protein
        return copy_protein

    def metropolis_accept(self, new_protein):
        """
        Calculate probability of accepting or rejecting move based on Metropolis criterion.
        --------
        Params:
            - new_protein (Protein): copied protein containing the insertion that has yet to be evaluated
        Returns:
            - Energy after insertion (avoids the need to calculate it again later on)
            - Probability of accepting
        """
        # Compute change of energy
        E_before = self.compute_energy(self.protein.pose)
        E_after = self.compute_energy(new_protein.pose)

        E_change = E_after - E_before

        # Obtain probability based on Metropolis criteria
        if E_change <= 0:
            return [E_after, 1]
        else:
            return [E_after, np.exp(-E_change / self.curr_T)]


    def anneal_temp(self):
        """
        Anneal temperature using exponential annealing schedule. Consider kT to be a single variable (i.e. ignore Boltzmann constant)
        --------
        Updates the temperature after step is accepted, otherwise it is not called
        (e.g. when fragment is rejected or fragment sample space is exhausted and new position is sampled)
        """
        self.curr_T = self.annealing_rate*self.curr_T

    def step(self,position = None):
        """
        Take a single MCMC step. Each step should do the following:
        1. sample position in chain
            - Note: think about positions you can sample a k-mer fragment from. 
        2. sample fragment at that position and replace torsions in a *copied version* of the protein
        3. measure energy after replacing fragment
        4. accept or reject based on Metropolis criterion
            - if accept: incorporate proposed insertion and anneal temperature
            - if reject: sample new fragment (go to step 2)
        """

        # Checks if we are recursively calling function with a specified position
        # or want to sample a new position
        if position == None:

            random.seed()
            position = random.randint(self.pos_range[0],self.pos_range[1])

        # Copy protein and insert fragment from sample spacce
        new_protein = self.perturb_fragment(position)

        # Evaluate insertion by Metropolis critera
        [E_after, probability] = self.metropolis_accept(new_protein)

        # Generate new seed and sample uniformly from [0,1]
        random.seed()
        value = random.uniform(0, 1)

        # If uniform sampling <= probability, accept
        if value <= probability:

            # Update protein
            self.protein = new_protein

            # Update temperature
            self.anneal_temp()

            # Reinitialize set containing fragments sampled at this step
            self.sampled_candidates = set()

            # Return energy of updated protein
            return E_after

        # If we exahused fragment sample space without accepting
        elif self.sampled_candidates == set(range(self.N)):

            # Reinitialize set containing fragments sampled at this step
            self.sampled_candidates = set()

            return None

        # Otherwise, recursively call to sample new fragments at this position until
        # fragment sample space is exhausted (second case) or fragment is accepted (first case)
        else:
            return self.step(position)


    def simulate(self):
        """
        Run full MCMC simulation from start_temp to end_temp. 
        Be sure to save the best (lowest-energy) structure, so you can access it after.
        It is also a good idea to track certain variables during the simulation (temp, energy, and more).
        -------- 
        Returns:
            - Protein with the lowest energy configuration
        """

        # Initialize best protein and minimum energy
        best_protein = Protein(pose = self.protein.pose)
        min_energy = self.compute_energy(best_protein.pose)

        # Keep track of iterations
        iterations = 0

        # Keep track of data at each iteration; will use to create log file
        data = []
        
        # Iterate until temperature at step is <= end temperature
        while self.curr_T > self.T_end:

            # Get output (either protein (if insertion was accepted) or None)
            output = self.step()

            # If all fragments from sample space were rejected
            # do another step at a random position
            if output == None:
                continue
            else:

                # Append iteration data
                iterations += 1
                iteration_data = [str(iterations), str(self.curr_T), str(output)]
                data.append(iteration_data)

                # If the update protein had a lower energy, save it
                if output <= min_energy:
                    min_energy = output
                    best_protein = self.protein


        # Log iteration, temperature, and energy @ each iteration in the appropriately numbered simulation folder
        filename = f'{self.sim_folder}/log_{self.frag_length}mer'  

        with open(filename, 'w') as f:
            for entry in data:
                    f.write('\t'.join(entry[:]) + '\n')

        # Save the lowest-energy (i.e., best) protein after 3-mer refinement
        if self.frag_length == 3:
            best_protein.save_pdb(f'{self.sim_folder}/best.pdb')

        # Return best protein
        return best_protein







    



