from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB import *
import numpy as np
import random
from tqdm import tqdm

from pyrosetta import *
from pyrosetta.rosetta.core.scoring.constraints import AtomPairConstraint , DihedralConstraint

from pyrosetta.rosetta.core.kinematics import MoveMap
# from pyrosetta.rosetta.protocols.minimization_packing import MinMover

from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.core.scoring.func import CircularHarmonicFunc ,HarmonicFunc

from pyrosetta.rosetta.core.id import AtomID

from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.moves import MonteCarlo


from pyrosetta.rosetta.protocols.moves import RepeatMover



NUM_ITERATIONS = 10
STD_ANGLE  = 5
STD_BOUND = 0.01

# Define hyperparameters
monte_carlo_temp = 2.0  # Increase if acceptance rate is too low
num_relax_rounds = 2  # Increase for more thorough, but slower, optimization
max_relax_iter = 200  # Increase for more thorough, but slower, optimization
num_iterations = 100  # Increase for longer, but slower, runs


# Load structure
parser = PDBParser()
structure = parser.get_structure('Your_protein', 'test.pdb')

# Get chain A
chain_A = structure[0]['A']



# Save chain A as a separate file
io = PDBIO()
io.set_structure(chain_A)
io.save('base_line_chain_A.pdb')



# Extract sequence
ppb = PPBuilder()
sequence = ""
for pp in ppb.build_peptides(chain_A):
    sequence += str(pp.get_sequence())

print("Orginal Chain Sequence" , sequence)


# Filter out only alpha carbon atoms
alpha_carbons = [atom for residue in chain_A for atom in residue if atom.get_id() == "CA"]

# Define function to calculate distance
def calculate_distance(atom1, atom2):
    diff_vector  = atom2.coord - atom1.coord
    return np.sqrt(np.sum(diff_vector * diff_vector))

# Create distance matrix
num_atoms = len(alpha_carbons)
distance_matrix = np.zeros((num_atoms, num_atoms))

for i in range(num_atoms):
    for j in range(num_atoms):
        distance_matrix[i, j] = calculate_distance(alpha_carbons[i], alpha_carbons[j])

# Save distance matrix to a csv file
np.savetxt("distance_matrix.csv", distance_matrix, delimiter=",")


print("Shape Distance Matrix :" , distance_matrix.shape )



ppb = PPBuilder()

# Using the method you mentioned to get the list of phi and psi angles
phi_psi_list = ppb.build_peptides(chain_A)[0].get_phi_psi_list()

# Separate the phi and psi angles into two lists
phi_angles = [angles[0] for angles in phi_psi_list if angles[0] is not None]
psi_angles = [angles[1] for angles in phi_psi_list if angles[1] is not None]
import math

phi_angles_deg = [math.degrees(angle) for angle in phi_angles]
psi_angles_deg = [math.degrees(angle) for angle in psi_angles]

phi_stdevs = [random.uniform(0, STD_ANGLE) for _ in range(len(phi_angles_deg))]
psi_stdevs = [random.uniform(0, STD_ANGLE) for _ in range(len(psi_angles_deg))]

print(len(phi_angles_deg) , len(psi_angles_deg) , max(phi_angles_deg))


# initialize PyRosetta
init()

# load your protein
pose = pose_from_sequence(sequence)

# Assume phi_angles and psi_angles are your precalculated lists of angles, which are of length N-1

# Set the phi angles for residues 2 to N

for i in range(2, pose.total_residue() + 1):
    pose.set_phi(i, phi_angles_deg[i-2])

# Set the psi angles for residues 1 to N-1
for i in range(1, pose.total_residue()):
    pose.set_psi(i, psi_angles_deg[i-1])

for i in range(len(distance_matrix)):
    for j in range(i+1, len(distance_matrix[i])):
        # create an AtomPairConstraint for the alpha carbons of residues i+1 and j+1
        atom1 = AtomID(2, i+1)  # 2 is the index for alpha carbon ('CA') in Rosetta atom numbering
        atom2 = AtomID(2, j+1)
        func = HarmonicFunc(distance_matrix[i][j], STD_BOUND)  # Harmonic potential with standard deviation of 1.0
        constraint = AtomPairConstraint(atom1, atom2, func)

        # add the constraint to the pose
        pose.add_constraint(constraint)
  
      
for i in range(1, pose.total_residue()):  # Rosetta uses 1-indexing
    # Create AtomID for the atoms forming the psi angle: N(i)-CA(i)-C(i)-N(i+1)
    atom1 = AtomID(1, i)   # N atom of current residue
    atom2 = AtomID(2, i)   # CA atom of current residue
    atom3 = AtomID(4, i)   # C atom of current residue
    atom4 = AtomID(1, i+1) # N atom of next residue if it exists
    
    # Create a CircularHarmonicFunc for the psi angle
    func = CircularHarmonicFunc(psi_angles_deg[i-1], psi_stdevs[i-1])

    # Create a DihedralConstraint for the psi angle
    constraint = DihedralConstraint(atom1, atom2, atom3, atom4, func)

    # Add the constraints to the pose
    pose.add_constraint(constraint)

for i in range(2, pose.total_residue() + 1):  # Rosetta uses 1-indexing
    # Create AtomID for the atoms forming the phi angle: C(i-1)-N(i)-CA(i)-C(i)
    atom1 = AtomID(4, i-1)  # C atom of previous residue
    atom2 = AtomID(1, i)  # N atom of current residue
    atom3 = AtomID(2, i)  # CA atom of current residue
    atom4 = AtomID(4, i)  # C atom of current residue
    
    # Create a CircularHarmonicFunc for the phi angle
    func = CircularHarmonicFunc(phi_angles_deg[i-2], phi_stdevs[i-2])

    # Create a DihedralConstraint for the phi angle
    constraint = DihedralConstraint(atom1, atom2, atom3, atom4, func )

    # Add the constraints to the pose
    pose.add_constraint(constraint)
    

# exit()
# Create a MoveMap that will allow backbone torsion angles to change
movemap = MoveMap()
movemap.set_bb(True)  # True means all backbone torsion angles are allowed to change
movemap.set_chi(True)
# Create a score function
scorefxn = get_score_function()  # Use the default score function
from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep, fa_sol, fa_intra_rep, fa_elec, pro_close, hbond_sr_bb, hbond_lr_bb, hbond_bb_sc, hbond_sc, constraints

# Increase the weights of the constraint terms
scorefxn.set_weight(ScoreType.atom_pair_constraint, 10.0)
scorefxn.set_weight(ScoreType.dihedral_constraint, 10.0)

#consider other functions

scorefxn.set_weight(fa_atr, 0.8)
scorefxn.set_weight(fa_rep, 0.44)
scorefxn.set_weight(fa_sol, 0.75)
scorefxn.set_weight(fa_intra_rep, 0.004)
scorefxn.set_weight(fa_elec, 0.70)
scorefxn.set_weight(pro_close, 1.25)
scorefxn.set_weight(hbond_sr_bb, 0.5)
scorefxn.set_weight(hbond_lr_bb, 0.5)
scorefxn.set_weight(hbond_bb_sc, 0.5)
scorefxn.set_weight(hbond_sc, 1.0)



# Create a Monte Carlo mover
mc = MonteCarlo(pose, scorefxn, monte_carlo_temp)  # The last parameter is the temperature

# Create a FastRelax mover
fast_relax = FastRelax(scorefxn, num_relax_rounds)  # The second parameter is the number of rounds
fast_relax.set_movemap(movemap)
fast_relax.max_iter(max_relax_iter)
# fast_relax.max_iter(100000000)  # Set maximum iterations
fast_relax.set_scorefxn(scorefxn)


for _ in tqdm(range(NUM_ITERATIONS)):  # The number of iterations
    fast_relax.apply(pose)
    mc.boltzmann(pose)


pose.dump_pdb("output_chain_a.pdb")