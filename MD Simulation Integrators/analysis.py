from pymol import cmd
import matplotlib.pyplot as plt
import numpy as np

def get_residue_numbers(sel):
    """
    Return a list of integers giving residue numbers for all protein
    residues in the given selection.
    """
    resnums = {'resnums': []}
    cmd.iterate(f'{sel} and name CA', 'resnums.append(resi)', space=resnums)
    resnums = resnums['resnums']
    return [int(resnum) for resnum in resnums]

def get_n_states(sel):
    """
    Returns the number of states for the given selection.
    """
    return cmd.count_states(sel)

def plot_distance(entry, sel1, sel2):
    n_states = get_n_states(entry)
    cmd.distance(f'{entry} and {sel1}', f'{entry} and {sel2}')
    distances = []
    for state in range(1, n_states+1):
        distances += [cmd.get_distance(f'{entry} and {sel1}', f'{entry} and {sel2}', state=state)]

    plt.plot(distances)
    plt.ylabel(f'Distance from {sel1} to {sel2} (Angstroms)')
    plt.xlabel('Time')
    plt.title(entry)
    plt.show()

def plot_rmsd(ref_sel, query_sel1, query_sel2):
    n_states1 = get_n_states(query_sel1)
    rmsds1 = []
    for	state in range(1, n_states1+1):
        rmsds1 += [cmd.align(query_sel1, ref_sel, target_state=0, mobile_state=state, cycles=0, transform=0)[0]]

    t = np.linspace(0, 1000000, num=len(rmsds1))
    plt.plot(t,rmsds1, label="Reverse Simulation")
    plt.ylabel(f'RMSD (Angstroms)')
    plt.xlabel('Time')


    n_states2 = get_n_states(query_sel2)
    rmsds2 = []
    for state in range(1, n_states2+1):
        rmsds2 += [cmd.align(query_sel2, ref_sel, target_state=0, mobile_state=state, cycles=0, transform=0)[0]]

    plt.plot(t,rmsds2, label="Forward Simulation")
    plt.ylabel(f'RMSD (Angstroms)')
    plt.xlabel('Time')

    plt.legend()
    plt.show()

def ramachandran(sel):
    """
    Produce a Ramachandran plot for residues in the given selection.
    """
    cmd.delete('phi')
    cmd.delete('psi')
    resnums = get_residue_numbers(sel)
    phis, psis = [], []

    # For each residue in resnums, add the phi and psi angle to phis and psis, respectively.

    # PyMol has two commands related to dihedral angles:
    # cmd.dihedral(name, sel1, sel2, sel3, sel4) will plot the dihedral on the
    #     structure however it will not return the value.
    # cmd.get_dihedral(sel1, sel2, sel3, sel4) will return the dihedral but will
    #     not show it on the structure.

    # Only cmd.get_dihedral is strictly required in your implementation, but
    # we highly recommend that you call both commands with the same selections
    # so that you can visually see the angles for debugging purposes.
    # Note that cmd.dihedral has an additional ``name'' argument, which you
    # should set to "phi" for the phi angles and psi for the "psi" angles.

    # Some tips on what various error messages mean:
    # "Error: Selection 1: Not found": The first selection matches no atoms.
    # "Error: Selection 1: Invalid selection name": The first selection matches multiple atoms.
    # Equivalent messages for Selection 2 mean the second selection is invalid, and so on.

    ############################################################################
    # Edit here.

    # Reinitialize to avoid selection error
    cmd.reinitialize()

    # Parse filename and residue range from cmd input
    cmd_input = sel.split(" ")
    file = cmd_input[0]

    # Default selection
    selection = "all"

    # If selection is specified in input, update 'selection'
    if len(cmd_input) > 2:
        selection = " ".join(cmd_input[2:])

    # Load file
    cmd.load(f"pdbs/{file}.pdb")

    # Get list of residues
    residues = get_residue_numbers(selection)

    # Iterate through residues and compute phi, psi
    for i in range(0,len(residues)):

        if i + 2 < len(residues):

            ### USE TO ANSWER QUESTION #6 by reporting observations
            # Calculate omega (between Cαi − Ci − Ni+1 − Cαi+1)
            # print(cmd.get_dihedral(f"{residues[i+1]}/ca",f"{residues[i+1]}/c", f"{residues[i+2]}/n", f"{residues[i+2]}/ca"))

            # Calculate psi (between Ni − Cαi − Ci − Ni+1): get_dihedral (i)/n, (i)/ca, (i)/c, (i+1)/n
            psi = cmd.get_dihedral(f"{residues[i+1]}/n",f"{residues[i+1]}/ca", f"{residues[i+1]}/c", f"{residues[i+2]}/n")
            psis.append(psi)

            # Here, (i+1) if statement is contained in (i+2) because the pairs are generated together
            if i + 1 < len(residues):
                # Calculate phi (between Ci-1 − Ni − Cαi − Ci): get_dihedral (i-1)/c, (i)/n, (i)/ca, (i)/c : 
                phi = cmd.get_dihedral(f"{residues[i]}/c",f"{residues[i+1]}/n", f"{residues[i+1]}/ca", f"{residues[i+1]}/c")
                phis.append(phi)

    ############################################################################
    plt.scatter(phis, psis)
    plt.xlabel('phi')
    plt.ylabel('psi')
    plt.ylim(-180, 180)
    plt.xlim(-180, 180)
    plt.gca().set_aspect('equal')
    plt.show()

cmd.extend('plot_distance', plot_distance)
cmd.extend('plot_rmsd', plot_rmsd)
cmd.extend('ramachandran', ramachandran)
