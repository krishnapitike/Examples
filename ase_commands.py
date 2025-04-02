#################################################################
# reading and writing lammps data
#################################################################

from ase.io.lammpsdata import read_lammps_data, write_lammps_data
from ase.build import make_supercell

# reading lammps data
unit_cell=read_lammps_data('init.data',style='atomic')

# making supercell
transformation = [[2,0,0],[0,2,0],[0,0,2]]
super_cell=make_supercell(unit_cell, transformation)

# writing lammps data
write_lammps_data('supercell.data',super_cell,atom_style='atomic')

#################################################################
# sorting in ase based on elements
#################################################################

from ase.io import read, write
from ase import Atoms
import numpy as np
import random

# Desired element order
desired_order = ['Fe', 'Ni', 'Cr', 'Al']

# Load POSCAR
atoms = read('POSCAR')

# Find indices of Fe atoms
fe_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'Fe']

# Set the number of substitutions (e.g., 10% of Fe atoms)
num_substitute = int(0.1 * len(fe_indices))

# Choose random indices to substitute
substitute_indices = random.sample(fe_indices, num_substitute)

# Substitute Fe with Ni
for idx in substitute_indices:
    atoms[idx].symbol = 'Ni'

# Sort atoms by the desired element order
atoms_sorted = Atoms([atom for symbol in desired_order for atom in atoms if atom.symbol == symbol],
                     cell=atoms.cell, pbc=atoms.pbc)

# Write the modified structure to a new POSCAR
write('POSCAR_modified', atoms_sorted, format='vasp')
