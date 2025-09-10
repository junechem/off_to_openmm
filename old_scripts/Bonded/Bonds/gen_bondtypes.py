#!/usr/bin/env python3

def read_atoms(atoms_string):
   atoms = {}
   for line in atoms_string.strip().split('\n'):
       parts = line.split()
       if len(parts) >= 2:
           atom_number = int(parts[0])
           atom_type = parts[1]
           atoms[atom_number] = atom_type
   return atoms

def read_bonds(bonds_string):
   bonds = []
   for line in bonds_string.strip().split('\n'):
       parts = line.split()
       if len(parts) >= 2:
           atom1 = int(parts[0])
           atom2 = int(parts[1])
           bonds.append((atom1, atom2))
   return bonds

def group_bonds_by_type(bonds, atoms):
   bond_types = {}

   for bond in bonds:
       atom1, atom2 = bond

       # Get atom types for the bond
       type1 = atoms[atom1]
       type2 = atoms[atom2]

       # Create a consistent representation of the bond type
       # Always put the alphabetically smaller type first
       if type1 <= type2:
           type_signature = f"{type1}-{type2}"
       else:
           type_signature = f"{type2}-{type1}"

       # Add to the corresponding group
       if type_signature not in bond_types:
           bond_types[type_signature] = []

       bond_types[type_signature].append(bond)

   return bond_types

def process_bond_data(atoms_string, bonds_string):
   # Parse input data
   atoms = read_atoms(atoms_string)
   bonds = read_bonds(bonds_string)

   # Group bonds by type
   bond_types = group_bonds_by_type(bonds, atoms)

   # Write results to file
   with open('bondtypes.dat', 'w') as f:
       for bond_type, bonds_list in sorted(bond_types.items()):
           f.write(f"\n;Bond Type: {bond_type}\n")
           num_bonds = len(bonds_list)
           f.write(f"\tHAR {num_bonds}\tFIT 1.0\t700\n")
           for bond in sorted(bonds_list):
               f.write(f"\t{bond[0]} {bond[1]}\n")

with open('atoms.dat', 'r') as f:
    atoms_data = f.read()

with open('bonds.dat', 'r') as f:
    bonds_data = f.read()


process_bond_data(atoms_data, bonds_data)
