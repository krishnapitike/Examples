from ase.io import read, write
from ase.io.lammpsdata import read_lammps_data, write_lammps_data
from ase.neighborlist import NeighborList
from ase import Atoms, Atom
import numpy as np
import json
import random

# === Load parameters ===
with open("input.json", "r") as f:
    params = json.load(f)

input_format = params["input_format"]
output_format = params["output_format"]
filename = params["input_file"]
void_file = params["void_file"]
bubble_file = params["bubble_file"]
center = np.array(params["center"])
radius = float(params["radius"])
hemisphere = params.get("hemisphere", "top").lower()  # "top", "bottom", "full"
desired_order = params.get("order", ["W", "Ni", "Fe", "He"])  # Order of elements in the output
he_vacancy_ratio = params.get("he_vacancy_ratio", 4.0)  # Default: 4 He per vacancy
tetra_dist = 1.0


# === Read structure ===
if input_format == 'vasp':
    atoms = read(filename)
elif input_format == 'lammps_data':
    atoms = read_lammps_data(filename, atom_style='atomic')
else:
    raise ValueError(f"Unsupported input format: {input_format}")

cell = atoms.get_cell()
positions = atoms.get_positions()

# === Add dummy atom at the center ===
dummy_atom = Atom('X', position=center)
atoms_with_dummy = atoms + Atoms([dummy_atom])
dummy_index = len(atoms_with_dummy) - 1

# === Build NeighborList object (just once) ===
cutoffs = [radius] * len(atoms_with_dummy)
nl = NeighborList(cutoffs=cutoffs, skin=0.0, self_interaction=False, bothways=True)
nl.update(atoms_with_dummy)

# === Query neighbors of the dummy atom only ===
indices, offsets = nl.get_neighbors(dummy_index)

to_delete = set()
for idx, offset in zip(indices, offsets):
    pos = atoms_with_dummy[idx].position + offset @ cell
    if np.linalg.norm(pos - center) <= radius:
        if hemisphere == "top" and pos[2] > center[2]:
            to_delete.add(idx)
        elif hemisphere == "bottom" and pos[2] < center[2]:
            to_delete.add(idx)
        elif hemisphere == "full":
            to_delete.add(idx)

# === Control number of He atoms via He/vacancy ratio ===
num_vacancies = len(to_delete)
max_he = 4 * num_vacancies

# Target number of He atoms
n_he_desired = int(he_vacancy_ratio * num_vacancies)
n_he_desired = min(n_he_desired, max_he)  # Clip to max available

# === Define tetrahedral He offsets ===
tetra_offsets = tetra_dist * np.array([
    [1, 1, 1],
    [-1, -1, 1],
    [-1, 1, -1],
    [1, -1, -1]
]) / np.sqrt(3) / 1.63299

# === Generate He atoms ===
he_positions = []
for idx in to_delete:
    atom_pos = atoms[idx].position
    for offset in tetra_offsets:
        he_pos = atom_pos + offset
        frac = np.linalg.solve(cell.T, he_pos) % 1.0
        wrapped = cell.T @ frac
        he_positions.append(wrapped)

if n_he_desired < len(he_positions):
    he_positions = random.sample(he_positions, n_he_desired)
    print(f"Reduced He atoms from {len(he_positions)} to {n_he_desired} to match He/vacancy ratio = {he_vacancy_ratio}")
elif he_vacancy_ratio > 4.0:
    print(f"Warning: he_vacancy_ratio > 4.0 requested, but only 4 He atoms per vacancy available.")

print(f"Deleted {len(to_delete)} atoms from the {hemisphere} hemisphere and added {len(he_positions)} He atoms.")

# === Remove deleted atoms and add He atoms ===
keep_indices = [i for i in range(len(atoms)) if i not in to_delete]
atoms_trimmed = atoms[keep_indices]

he_atoms = Atoms('He' * len(he_positions),
                 positions=he_positions,
                 cell=cell,
                 pbc=True)

final_atoms = atoms_trimmed + he_atoms

final_atoms_sorted = Atoms([atom for symbol in desired_order for atom in final_atoms if atom.symbol == symbol],
                     cell=atoms.cell, pbc=atoms.pbc)

#write(output_file+'.vasp', final_atoms_sorted, direct=1)
#final_atoms_sorted = read(output_file+'.vasp')

if output_format == 'lammps_data':
    write_lammps_data(bubble_file, final_atoms_sorted, specorder=desired_order, atom_style='atomic',masses=True)
    write_lammps_data(void_file, atoms_trimmed, specorder=desired_order[:-1], atom_style='atomic',masses=True) 
elif output_format == 'vasp':
    write(bubble_file+'.vasp', final_atoms_sorted, direct=1)
    write(void_file+'.vasp', atoms_trimmed, direct=1)
elif output_format == 'both':
    write_lammps_data(bubble_file+'.vasp', final_atoms_sorted, specorder=desired_order, atom_style='atomic',masses=True)
    write(bubble_file+'.vasp', final_atoms_sorted, direct=1)
    write_lammps_data(void_file, atoms_trimmed, specorder=desired_order[:-1], atom_style='atomic',masses=True)
    write(void_file+'.vasp', atoms_trimmed, direct=1)
else:
    raise ValueError(f"Invalid output format: {output_format}")
