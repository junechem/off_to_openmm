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

def find_angles(bonds):
    # Create a dictionary mapping each atom to its connected atoms
    connections = {}
    for atom1, atom2 in bonds:
        if atom1 not in connections:
            connections[atom1] = []
        if atom2 not in connections:
            connections[atom2] = []
        connections[atom1].append(atom2)
        connections[atom2].append(atom1)
    
    # Find all angles
    angles = []
    for central_atom in connections:
        if len(connections[central_atom]) >= 2:
            connected_atoms = connections[central_atom]
            for i in range(len(connected_atoms)):
                for j in range(i+1, len(connected_atoms)):
                    angle = (connected_atoms[i], central_atom, connected_atoms[j])
                    angles.append(angle)
    
    return angles

def group_angles_by_type(angles, atoms):
    angle_types = {}
    
    for angle in angles:
        atom1, central_atom, atom3 = angle
        
        # Get atom types for this angle
        type1 = atoms[atom1]
        type_central = atoms[central_atom]
        type3 = atoms[atom3]
        
        # Create a signature for this angle type
        type_signature = f"{type1}-{type_central}-{type3}"
        
        # Add to the corresponding group
        if type_signature not in angle_types:
            angle_types[type_signature] = []
        
        angle_types[type_signature].append(angle)
    
    return angle_types

def process_data(atoms_string, bonds_string):
    # Parse input data
    atoms = read_atoms(atoms_string)
    bonds = read_bonds(bonds_string)
    
    # Find all angles
    angles = find_angles(bonds)
    
    # Group angles by type
    angle_types = group_angles_by_type(angles, atoms)
    
    # Write results to file
    with open('angletypes.dat', 'w') as f:
        for angle_type, angles_list in angle_types.items():
            f.write(f"\n;Angle Type: {angle_type}\n")
            num_angles = len(angles_list)
            f.write(f"\t{'HAR':>4}{num_angles:>4}{'FIT':>4}{100:>4}{76:>4}\n")
            for angle in sorted(angles_list):
                f.write(f"\t\t{angle[0]} {angle[1]} {angle[2]}\n")

# Example usage with the provided data

with open('atoms.dat', 'r') as f:
    atoms_data = f.read()

with open('bonds.dat', 'r') as f:
    bonds_data = f.read()

process_data(atoms_data, bonds_data)
