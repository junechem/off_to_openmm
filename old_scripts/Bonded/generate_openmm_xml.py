#!/usr/bin/env python3

import math
import os

def parse_parameters(filename):
    """Parse the parameters.dat file and extract bond, angle, and dihedral parameters"""
    parameters = {
        'bonds': {},
        'angles': {},
        'dihedrals': {}
    }
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#define'):
                parts = line.split()
                if len(parts) >= 4:
                    param_name = parts[1]
                    param_type = parts[2]  # HAR or NCO
                    param_values = parts[3:]
                    
                    # Count underscores to determine type
                    underscore_count = param_name.count('_')
                    
                    if underscore_count == 1:  # Bond
                        parameters['bonds'][param_name] = {
                            'type': param_type,
                            'values': param_values
                        }
                    elif underscore_count == 2:  # Angle
                        parameters['angles'][param_name] = {
                            'type': param_type,
                            'values': param_values
                        }
                    elif underscore_count == 3:  # Dihedral
                        parameters['dihedrals'][param_name] = {
                            'type': param_type,
                            'values': param_values
                        }
    
    return parameters

def read_atoms(filename):
    """Read atoms.dat file and return atom mapping"""
    atoms = {}
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                atom_id = int(parts[0])
                atom_type = parts[1]
                atoms[atom_id] = atom_type
    return atoms

def read_unique_atom_types(filename):
    """Read atoms.dat file and return unique atom types"""
    unique_types = []
    
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                atom_type = parts[1]
                if atom_type not in unique_types:
                    unique_types.append(atom_type)
    
    return unique_types

def read_charges(filename):
    """Read charges.dat file and return charge mapping by atom type"""
    charges = {}
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                atom_type = parts[0]
                charge = float(parts[1])
                charges[atom_type] = charge
    return charges

def get_element_from_atom_type(atom_type):
    """Extract element from atom type (first letter)"""
    return atom_type[0]

def generate_atomtypes_xml(atom_types):
    """Generate AtomTypes XML section"""
    xml_lines = []
    xml_lines.append('<AtomTypes>')
    
    for atom_type in atom_types:
        element = get_element_from_atom_type(atom_type)
        xml_lines.append(f'<Type name="{atom_type}" class="{atom_type}" element="{element}" mass="0.0"/>')
    
    xml_lines.append('</AtomTypes>')
    
    return '\n'.join(xml_lines)

def create_unique_atom_names(atoms):
    """Create unique atom names based on element with sequential numbering"""
    atom_id_to_name = {}  # atom_id -> unique_name
    atom_id_to_type = {}  # atom_id -> atom_type
    element_counters = {}  # element -> count
    
    # Sort atom IDs to ensure consistent ordering
    for atom_id in sorted(atoms.keys()):
        atom_type = atoms[atom_id]
        element = atom_type[0]  # First letter is the element
        
        # Initialize counter for this element if not seen before
        if element not in element_counters:
            element_counters[element] = 0
        
        # Create unique name
        unique_name = f"{element}{element_counters[element]}"
        
        # Store mappings
        atom_id_to_name[atom_id] = unique_name
        atom_id_to_type[atom_id] = atom_type
        
        # Increment counter for this element
        element_counters[element] += 1
    
    return atom_id_to_name, atom_id_to_type

def read_bonds_from_bondtypes(filename):
    """Read actual bond pairs from bondtypes.dat file"""
    bonds = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith(';') and not line.startswith('HAR'):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        atom1_id = int(parts[0])
                        atom2_id = int(parts[1])
                        bonds.append((atom1_id, atom2_id))
                    except ValueError:
                        # Skip lines that don't have integer atom IDs
                        continue
    
    return bonds

def generate_residue_xml(atom_id_to_name, atom_id_to_type, bonds, charges):
    """Generate Residue XML section with unique atom names, explicit bonds, and charges"""
    xml_lines = []
    xml_lines.append('<Residues>')
    xml_lines.append('<Residue name="UNK">')
    
    # Generate <Atom> entries with unique names and charges
    for atom_id in sorted(atom_id_to_name.keys()):
        unique_name = atom_id_to_name[atom_id]
        atom_type = atom_id_to_type[atom_id]
        charge = charges.get(atom_type, 0.0)
        xml_lines.append(f'<Atom name="{unique_name}" type="{atom_type}" charge="{charge}"/>')
    
    # Generate <Bond> entries with unique atom names
    for atom1_id, atom2_id in bonds:
        if atom1_id in atom_id_to_name and atom2_id in atom_id_to_name:
            atom1_name = atom_id_to_name[atom1_id]
            atom2_name = atom_id_to_name[atom2_id]
            xml_lines.append(f'<Bond atomName1="{atom1_name}" atomName2="{atom2_name}"/>')
    
    xml_lines.append('</Residue>')
    xml_lines.append('</Residues>')
    
    return '\n'.join(xml_lines)

def process_bonds(bondtypes_file, parameters):
    """Process bonds from bondtypes.dat file"""
    bond_types = {}
    
    with open(bondtypes_file, 'r') as f:
        current_bond_type = None
        
        for line in f:
            line = line.strip()
            if line.startswith(';Bond Type:'):
                current_bond_type = line.split(':', 1)[1].strip()
            elif line and not line.startswith(';') and not line.startswith('HAR'):
                if current_bond_type:
                    # Extract atom types from bond_type (e.g., "C1-C2" -> ["C1", "C2"])
                    atom_types = current_bond_type.split('-')
                    if len(atom_types) == 2:
                        type1, type2 = atom_types
                        # Create normalized key for parameter lookup
                        param_key = f"{type1}_{type2}"
                        if param_key not in parameters['bonds']:
                            param_key = f"{type2}_{type1}"
                        
                        # Store bond type info
                        sorted_key = tuple(sorted([type1, type2]))
                        if sorted_key not in bond_types:
                            bond_types[sorted_key] = {
                                'param_key': param_key,
                                'type1': type1,
                                'type2': type2
                            }
    
    return bond_types

def process_angles(angletypes_file, parameters):
    """Process angles from angletypes.dat file"""
    angle_types = {}
    
    with open(angletypes_file, 'r') as f:
        current_angle_type = None
        
        for line in f:
            line = line.strip()
            if line.startswith(';Angle Type:'):
                current_angle_type = line.split(':', 1)[1].strip()
            elif line and not line.startswith(';') and not line.startswith('HAR'):
                if current_angle_type:
                    # Extract atom types from angle_type (e.g., "C2-C1-H1" -> ["C2", "C1", "H1"])
                    atom_types = current_angle_type.split('-')
                    if len(atom_types) == 3:
                        type1, type2, type3 = atom_types
                        # Create normalized key for parameter lookup
                        param_key = f"{type1}_{type2}_{type3}"
                        if param_key not in parameters['angles']:
                            param_key = f"{type3}_{type2}_{type1}"
                        
                        # Store angle type info
                        key = f"{type1}-{type2}-{type3}"
                        if key not in angle_types:
                            angle_types[key] = {
                                'param_key': param_key,
                                'type1': type1,
                                'type2': type2,
                                'type3': type3
                            }
    
    return angle_types

def process_dihedrals(dihedraltypes_file, parameters):
    """Process dihedrals from dihedraltypes.dat file"""
    dihedral_types = {}
    
    with open(dihedraltypes_file, 'r') as f:
        current_dihedral_type = None
        
        for line in f:
            line = line.strip()
            if line.startswith(';Dihedral Type:'):
                current_dihedral_type = line.split(':', 1)[1].strip()
            elif line and not line.startswith(';') and not line.startswith('NCO'):
                if current_dihedral_type:
                    # Extract atom types from dihedral_type (e.g., "C2-C2-C1-O0" -> ["C2", "C2", "C1", "O0"])
                    atom_types = current_dihedral_type.split('-')
                    if len(atom_types) == 4:
                        type1, type2, type3, type4 = atom_types
                        # Create normalized key for parameter lookup
                        param_key = f"{type1}_{type2}_{type3}_{type4}"
                        if param_key not in parameters['dihedrals']:
                            param_key = f"{type4}_{type3}_{type2}_{type1}"
                        
                        # Store dihedral type info
                        key = f"{type1}-{type2}-{type3}-{type4}"
                        if key not in dihedral_types:
                            dihedral_types[key] = {
                                'param_key': param_key,
                                'type1': type1,
                                'type2': type2,
                                'type3': type3,
                                'type4': type4
                            }
    
    return dihedral_types

def generate_xml(parameters, atom_types, atom_id_to_name, atom_id_to_type, bonds, bond_types, angle_types, dihedral_types, charges):
    """Generate the OpenMM XML format"""
    xml_content = []
    
    # Generate AtomTypes section (first)
    atomtypes_xml = generate_atomtypes_xml(atom_types)
    xml_content.append(atomtypes_xml)
    
    # Generate Residues section with unique atom names and charges
    residue_xml = generate_residue_xml(atom_id_to_name, atom_id_to_type, bonds, charges)
    xml_content.append(residue_xml)
    
    # Generate HarmonicBondForce section
    xml_content.append('<HarmonicBondForce>')
    for bond_key, bond_info in bond_types.items():
        param_key = bond_info['param_key']
        if param_key in parameters['bonds']:
            param_data = parameters['bonds'][param_key]
            if len(param_data['values']) >= 2:
                length = param_data['values'][0]
                k = param_data['values'][1]
                xml_content.append(f'<Bond type1="{bond_info["type1"]}" type2="{bond_info["type2"]}" length="{length}" k="{k}"/>')
    xml_content.append('</HarmonicBondForce>')
    
    # Generate HarmonicAngleForce section
    xml_content.append('<HarmonicAngleForce>')
    for angle_key, angle_info in angle_types.items():
        param_key = angle_info['param_key']
        if param_key in parameters['angles']:
            param_data = parameters['angles'][param_key]
            if len(param_data['values']) >= 2:
                theta_deg = param_data['values'][0]  # In degrees from parameters
                k = param_data['values'][1]
                
                # Convert angle from degrees to radians (OpenMM uses radians)
                try:
                    theta_deg_float = float(theta_deg)
                    theta_rad = theta_deg_float * math.pi / 180.0
                    xml_content.append(f'<Angle class1="{angle_info["type1"]}" class2="{angle_info["type2"]}" class3="{angle_info["type3"]}" angle="{theta_rad}" k="{k}"/>')
                except ValueError:
                    # If theta is not a number, use it as-is
                    xml_content.append(f'<Angle class1="{angle_info["type1"]}" class2="{angle_info["type2"]}" class3="{angle_info["type3"]}" angle="{theta_deg}" k="{k}"/>')
    xml_content.append('</HarmonicAngleForce>')
    
    # Generate PeriodicTorsionForce section
    xml_content.append('<PeriodicTorsionForce>')
    for dihedral_key, dihedral_info in dihedral_types.items():
        param_key = dihedral_info['param_key']
        if param_key in parameters['dihedrals']:
            param_data = parameters['dihedrals'][param_key]
            if len(param_data['values']) >= 3:
                phase = param_data['values'][0]  # Already in degrees, convert to radians
                k = param_data['values'][1]
                periodicity = param_data['values'][2]
                
                # Convert phase from degrees to radians (OpenMM uses radians)
                try:
                    phase_deg = float(phase)
                    phase_rad = phase_deg * math.pi / 180.0
                    xml_content.append(f'<Proper class1="{dihedral_info["type1"]}" class2="{dihedral_info["type2"]}" class3="{dihedral_info["type3"]}" class4="{dihedral_info["type4"]}" periodicity1="{periodicity}" phase1="{phase_rad}" k1="{k}"/>')
                except ValueError:
                    # If phase is not a number, use it as-is (might be 0.0000000)
                    xml_content.append(f'<Proper class1="{dihedral_info["type1"]}" class2="{dihedral_info["type2"]}" class3="{dihedral_info["type3"]}" class4="{dihedral_info["type4"]}" periodicity1="{periodicity}" phase1="{phase}" k1="{k}"/>')
    xml_content.append('</PeriodicTorsionForce>')
    
    return '\n'.join(xml_content)

def main():
    """Main function to generate OpenMM XML file"""
    
    # Parse parameters
    parameters = parse_parameters('parameters.dat')
    
    # Read atoms (using bonds directory as primary source)
    atoms = read_atoms('Bonds/atoms.dat')
    
    # Get unique atom types for AtomTypes section
    atom_types = read_unique_atom_types('Bonds/atoms.dat')
    
    # Read charges from charges.dat
    charges = read_charges('charges.dat')
    
    # Create unique atom names for Residue section
    atom_id_to_name, atom_id_to_type = create_unique_atom_names(atoms)
    
    # Read actual bonds from bondtypes.dat
    bonds = read_bonds_from_bondtypes('Bonds/bondtypes.dat')
    
    # Process bond types
    bond_types = process_bonds('Bonds/bondtypes.dat', parameters)
    
    # Process angle types
    angle_types = process_angles('Angles/angletypes.dat', parameters)
    
    # Process dihedral types
    dihedral_types = process_dihedrals('Dihedrals/dihedraltypes.dat', parameters)
    
    # Generate XML content
    xml_content = generate_xml(parameters, atom_types, atom_id_to_name, atom_id_to_type, bonds, bond_types, angle_types, dihedral_types, charges)
    
    # Write output file
    with open('forcefield.xml', 'w') as f:
        f.write(xml_content + '\n')

if __name__ == "__main__":
    main()