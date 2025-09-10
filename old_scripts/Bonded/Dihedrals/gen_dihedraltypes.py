#!/usr/bin/env python3
import os
import sys
### Usage: dihedrals.py (no arguments needed)
### Script looks for atoms.dat, bonds.dat, and valid_dihedrals.dat in current directory
### Atoms.dat contains atom numbers and atom types in rows
### Bonds.dat contains bonds (using atom indices) in rows
### valid_dihedrals.dat contains dihedrals (of form 'C2-C3-C3-C2') which should be printed at the
###     end of this script. One per row

# Default filenames
ATOMS_FILE = 'atoms.dat'
BONDS_FILE = 'bonds.dat'
VALID_DIHEDRALS_FILE = 'valid_dihedrals.dat'






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
    
    return angles, connections

def find_dihedrals(angles, connections):
    dihedrals = []
    
    # For each angle, check if we can extend it to a dihedral
    for angle in angles:
        atom1, atom2, atom3 = angle
        
        # Try to extend from atom1
        for atom0 in connections[atom1]:
            if atom0 != atom2:  # Make sure it's not already in our angle
                dihedral = (atom0, atom1, atom2, atom3)
                dihedrals.append(dihedral)
        
        # Try to extend from atom3
        for atom4 in connections[atom3]:
            if atom4 != atom2:  # Make sure it's not already in our angle
                dihedral = (atom1, atom2, atom3, atom4)
                dihedrals.append(dihedral)
    
    # Remove duplicates (considering dihedrals can be read in both directions)
    unique_dihedrals = []
    seen = set()
    
    for dihedral in dihedrals:
        # Create a canonical representation for this dihedral
        canonical = tuple(sorted([dihedral, dihedral[::-1]]))
        
        if canonical not in seen:
            seen.add(canonical)
            unique_dihedrals.append(dihedral)
    
    return unique_dihedrals

def get_dihedral_type(dihedral, atoms):
    return "-".join([atoms[atom] for atom in dihedral])

def normalize_dihedral_type(dihedral_type):
    """Create a canonical representation of dihedral type that considers reversals"""
    parts = dihedral_type.split("-")
    forward = "-".join(parts)
    reverse = "-".join(parts[::-1])
    return min(forward, reverse)

def is_valid_dihedral_type(dihedral_type, valid_types):
    """Check if a dihedral type matches any of the valid types, accounting for reversals"""
    normalized_type = normalize_dihedral_type(dihedral_type)
    
    for valid_type in valid_types:
        normalized_valid = normalize_dihedral_type(valid_type)
        if normalized_type == normalized_valid:
            return True
    
    return False

def group_dihedrals_by_type(dihedrals, atoms, valid_dihedral_types=None):
    # First, normalize all valid dihedral types
    normalized_valid_types = set()
    if valid_dihedral_types:
        for vtype in valid_dihedral_types:
            normalized_valid_types.add(normalize_dihedral_type(vtype))
    
    dihedral_types = {}
    
    for dihedral in dihedrals:
        # Get atom types for this dihedral
        type_signature = get_dihedral_type(dihedral, atoms)
        normalized_signature = normalize_dihedral_type(type_signature)
        
        # Skip if we're filtering and this isn't a valid type
        if normalized_valid_types and normalized_signature not in normalized_valid_types:
            continue
        
        # Use the original signature for grouping (keeps the original order)
        if type_signature not in dihedral_types:
            dihedral_types[type_signature] = []
        
        dihedral_types[type_signature].append(dihedral)
    
    # Merge reversed types if both exist in the results
    merged_types = {}
    processed = set()
    
    for type_sig in dihedral_types:
        if type_sig in processed:
            continue
            
        normalized = normalize_dihedral_type(type_sig)
        reversed_sig = "-".join(type_sig.split("-")[::-1])
        
        if reversed_sig in dihedral_types and type_sig != reversed_sig:
            # Use the first one alphabetically as the canonical representation
            canonical = min(type_sig, reversed_sig)
            merged_types[canonical] = dihedral_types[type_sig] + dihedral_types[reversed_sig]
            processed.add(type_sig)
            processed.add(reversed_sig)
        else:
            merged_types[type_sig] = dihedral_types[type_sig]
            processed.add(type_sig)
    
    return merged_types

def process_data(atoms_string, bonds_string, valid_dihedral_types=None):
    # Parse input data
    atoms = read_atoms(atoms_string)
    bonds = read_bonds(bonds_string)
    
    # Find all angles and connections
    angles, connections = find_angles(bonds)
    
    # Find all possible dihedrals
    dihedrals = find_dihedrals(angles, connections)
    
    # Group dihedrals by type and filter if needed
    dihedral_types = group_dihedrals_by_type(dihedrals, atoms, valid_dihedral_types)
    
    # Write results to file
    with open('dihedraltypes.dat', 'w') as f:
        for dihedral_type, dihedrals_list in dihedral_types.items():
            f.write(f"\n;Dihedral Type: {dihedral_type}\n")
            num_dih = len(dihedrals_list)
            f.write(f"\t{'NCO':>4}{num_dih:>4}{'FIT':>4}{2.5:>4}{'X':>4}{'Y':>4}\n")
            for dihedral in sorted(dihedrals_list):
                f.write(f"\t\t{dihedral[0]} {dihedral[1]} {dihedral[2]} {dihedral[3]}\n")

def check_file_exists(filename):
    """Check if file exists and provide error message if not"""
    if not os.path.exists(filename):
        print(f"Error: Required file '{filename}' not found in current directory.")
        return False
    return True

def main():
    """Main function to run the dihedral analysis"""
    # Check if all required files exist
    required_files = [ATOMS_FILE, BONDS_FILE, VALID_DIHEDRALS_FILE]
    missing_files = []
    
    for filename in required_files:
        if not check_file_exists(filename):
            missing_files.append(filename)
    
    if missing_files:
        print("\nThe following files are required in the current directory:")
        print(f"- {ATOMS_FILE}: contains atom numbers and atom types in rows")
        print(f"- {BONDS_FILE}: contains bonds (using atom indices) in rows") 
        print(f"- {VALID_DIHEDRALS_FILE}: contains valid dihedrals (e.g., 'C2-C3-C3-C2') one per row")
        sys.exit(1)
    
    # Read the files
    try:
        with open(ATOMS_FILE, 'r') as f:
            atoms_data = f.read()
        
        with open(BONDS_FILE, 'r') as f:
            bonds_data = f.read()
        
        with open(VALID_DIHEDRALS_FILE, 'r') as f:
            valid_dihedral_types = f.read().strip().split('\n')
            # Remove empty lines
            valid_dihedral_types = [line for line in valid_dihedral_types if line.strip()]
        
        # Process the data
        process_data(atoms_data, bonds_data, valid_dihedral_types)
        
    except Exception as e:
        print(f"Error reading files: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
