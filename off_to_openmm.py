#!/usr/bin/env python3
"""
Convert CRYOFF .off force field files to OpenMM .xml format.
"""

import argparse
import sys
import re
import math
from collections import defaultdict


class OffFileParser:
    """Parser for CRYOFF .off files with support for multiple format variations."""
    
    def __init__(self, filename):
        self.filename = filename
        self.molecules = {}
        self.parameters = {}
        self.nonbonded_interactions = {'COU': [], 'EXP': [], 'SRD': []}
        
    def parse(self):
        """Parse the .off file and extract all relevant data."""
        with open(self.filename, 'r') as f:
            content = f.read()
        
        # Parse molecule definitions
        self._parse_molecules(content)
        
        # Parse parameter definitions (#define statements)
        self._parse_parameters(content)
        
        # Parse nonbonded interactions
        self._parse_nonbonded_interactions(content)
        
        return self.molecules, self.parameters, self.nonbonded_interactions
    
    def _parse_molecules(self, content):
        """Parse molecule definitions from MOLTYP/MOL sections."""
        # Handle both [MOLTYP] and [MOL] formats
        mol_pattern = r'\[\s*(MOLTYP|MOL)\s*\]\s+(\w+)'
        
        lines = content.split('\n')
        i = 0
        
        while i < len(lines):
            line = lines[i].strip()
            mol_match = re.search(mol_pattern, line)
            
            if mol_match:
                mol_name = mol_match.group(2)
                self.molecules[mol_name] = {
                    'atoms': {},
                    'bonds': [],
                    'angles': [],
                    'dihedrals': [],
                    'exclusions': []
                }
                
                i = self._parse_molecule_section(lines, i + 1, mol_name)
            else:
                i += 1
    
    def _parse_molecule_section(self, lines, start_idx, mol_name):
        """Parse a single molecule section."""
        i = start_idx
        
        while i < len(lines):
            line = lines[i].strip()
            
            # Check for next molecule or end of section
            if re.search(r'\[\s*(MOLTYP|MOL)\s*\]', line) or line.startswith('[ COU ]'):
                break
                
            # Parse atoms section
            if re.search(r'\[\s*ATOMS?\s*\]', line):
                i = self._parse_atoms_section(lines, i + 1, mol_name)
            # Parse bonds section  
            elif re.search(r'\[\s*BONDS?\s*\]', line):
                i = self._parse_bonds_section(lines, i + 1, mol_name)
            # Parse angles section
            elif re.search(r'\[\s*ANGLES?\s*\]', line):
                i = self._parse_angles_section(lines, i + 1, mol_name)
            # Parse dihedrals section
            elif re.search(r'\[\s*DIHEDRALS?\s*\]', line):
                i = self._parse_dihedrals_section(lines, i + 1, mol_name)
            # Parse exclusions section
            elif re.search(r'\[\s*EXC\s*\]', line):
                i = self._parse_exclusions_section(lines, i + 1, mol_name)
            else:
                i += 1
                
        return i
    
    def _parse_atoms_section(self, lines, start_idx, mol_name):
        """Parse atoms section."""
        i = start_idx
        # Skip the atom count line
        if i < len(lines) and lines[i].strip().isdigit():
            i += 1
            
        while i < len(lines):
            line = lines[i].strip()
            
            if line.startswith('[') or not line:
                break
                
            parts = line.split()
            if len(parts) >= 3:
                try:
                    # Handle special atom IDs like "4*"
                    atom_id_str = parts[0].rstrip('*')
                    atom_id = int(atom_id_str)
                    atom_name = parts[1]
                    atom_type = parts[2]
                    self.molecules[mol_name]['atoms'][atom_id] = {
                        'name': atom_name,
                        'type': atom_type,
                        'special': '*' in parts[0]  # Mark if it's a special atom
                    }
                except ValueError:
                    # Skip lines that can't be parsed as atoms
                    pass
            i += 1
            
        return i
    
    def _parse_bonds_section(self, lines, start_idx, mol_name):
        """Parse bonds section."""
        i = start_idx
        # Skip the bond count line
        if i < len(lines) and lines[i].strip().isdigit():
            i += 1
            
        current_bond_group = None
        
        while i < len(lines):
            line = lines[i].strip()
            
            if line.startswith('[') or not line:
                break
                
            if line.startswith('HAR'):
                # New bond group definition
                current_bond_group = line
                i += 1
                continue
                
            # Bond atom pairs
            parts = line.split()
            if len(parts) >= 2 and current_bond_group:
                try:
                    atom1_id = int(parts[0])
                    atom2_id = int(parts[1])
                    self.molecules[mol_name]['bonds'].append({
                        'atoms': (atom1_id, atom2_id),
                        'definition': current_bond_group
                    })
                except ValueError:
                    pass
            i += 1
            
        return i
    
    def _parse_angles_section(self, lines, start_idx, mol_name):
        """Parse angles section."""
        i = start_idx
        # Skip the angle count line
        if i < len(lines) and lines[i].strip().isdigit():
            i += 1
            
        current_angle_group = None
        
        while i < len(lines):
            line = lines[i].strip()
            
            if line.startswith('[') or not line:
                break
                
            if line.startswith('HAR'):
                # New angle group definition
                current_angle_group = line
                i += 1
                continue
                
            # Angle atom triplets
            parts = line.split()
            if len(parts) >= 3 and current_angle_group:
                try:
                    atom1_id = int(parts[0])
                    atom2_id = int(parts[1])
                    atom3_id = int(parts[2])
                    self.molecules[mol_name]['angles'].append({
                        'atoms': (atom1_id, atom2_id, atom3_id),
                        'definition': current_angle_group
                    })
                except ValueError:
                    pass
            i += 1
            
        return i
    
    def _parse_dihedrals_section(self, lines, start_idx, mol_name):
        """Parse dihedrals section."""
        i = start_idx
        # Skip the dihedral count line
        if i < len(lines) and lines[i].strip().isdigit():
            i += 1
            
        current_dihedral_group = None
        
        while i < len(lines):
            line = lines[i].strip()
            
            if line.startswith('[') or not line:
                break
                
            if line.startswith('NCO'):
                # New dihedral group definition
                current_dihedral_group = line
                i += 1
                continue
                
            # Dihedral atom quadruplets
            parts = line.split()
            if len(parts) >= 4 and current_dihedral_group:
                try:
                    atom1_id = int(parts[0])
                    atom2_id = int(parts[1])
                    atom3_id = int(parts[2])
                    atom4_id = int(parts[3])
                    self.molecules[mol_name]['dihedrals'].append({
                        'atoms': (atom1_id, atom2_id, atom3_id, atom4_id),
                        'definition': current_dihedral_group
                    })
                except ValueError:
                    pass
            i += 1
            
        return i
    
    def _parse_exclusions_section(self, lines, start_idx, mol_name):
        """Parse exclusions section."""
        i = start_idx
        # Skip the exclusion count line
        if i < len(lines) and lines[i].strip().isdigit():
            i += 1
            
        while i < len(lines):
            line = lines[i].strip()
            
            if line.startswith('[') or not line:
                break
                
            parts = line.split()
            if len(parts) >= 2:
                try:
                    atom1_id = int(parts[0])
                    atom2_id = int(parts[1])
                    self.molecules[mol_name]['exclusions'].append((atom1_id, atom2_id))
                except ValueError:
                    pass
            i += 1
            
        return i
    
    def _parse_parameters(self, content):
        """Parse #define parameter statements."""
        define_pattern = r'#define\s+(\S+)\s+(\S+)\s+(.*)'
        
        for line in content.split('\n'):
            match = re.search(define_pattern, line.strip())
            if match:
                param_name = match.group(1)
                param_type = match.group(2)
                param_values = match.group(3).split()
                
                # Categorize by number of underscores
                underscore_count = param_name.count('_')
                
                if underscore_count == 1:  # Bond
                    if 'bonds' not in self.parameters:
                        self.parameters['bonds'] = {}
                    self.parameters['bonds'][param_name] = {
                        'type': param_type,
                        'values': param_values
                    }
                elif underscore_count == 2:  # Angle
                    if 'angles' not in self.parameters:
                        self.parameters['angles'] = {}
                    self.parameters['angles'][param_name] = {
                        'type': param_type,
                        'values': param_values
                    }
                elif underscore_count == 3:  # Dihedral
                    if 'dihedrals' not in self.parameters:
                        self.parameters['dihedrals'] = {}
                    self.parameters['dihedrals'][param_name] = {
                        'type': param_type,
                        'values': param_values
                    }
    
    def _parse_nonbonded_interactions(self, content):
        """Parse COU, EXP, and SRD interaction sections."""
        lines = content.split('\n')
        i = 0
        
        while i < len(lines):
            line = lines[i].strip()
            
            # Parse COU section
            if line.startswith('[ COU ]'):
                i = self._parse_interaction_section(lines, i + 1, 'COU')
            # Parse EXP section
            elif line.startswith('[ EXP ]'):
                i = self._parse_interaction_section(lines, i + 1, 'EXP')
            # Parse SRD section
            elif line.startswith('[ SRD ]'):
                i = self._parse_interaction_section(lines, i + 1, 'SRD')
            else:
                i += 1
    
    def _parse_interaction_section(self, lines, start_idx, interaction_type):
        """Parse a nonbonded interaction section."""
        i = start_idx
        # Skip the interaction count line
        if i < len(lines) and lines[i].strip().isdigit():
            i += 1
            
        while i < len(lines):
            line = lines[i].strip()
            
            if line.startswith('[') or not line:
                break
                
            parts = line.split()
            if len(parts) >= 4:
                atom1 = parts[0]
                atom2 = parts[1]
                fix_flag = parts[2]
                values = parts[3:]
                
                self.nonbonded_interactions[interaction_type].append({
                    'atom1': atom1,
                    'atom2': atom2,
                    'fix_flag': fix_flag,
                    'values': values
                })
            i += 1
            
        return i


def read_charges_file(filename):
    """Read charges file with 'AtomType Charge' format."""
    charges = {}
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    atom_type = parts[0]
                    charge = float(parts[1])
                    charges[atom_type] = charge
    except FileNotFoundError:
        print(f"Warning: Charges file '{filename}' not found")
        return {}
    except Exception as e:
        print(f"Error reading charges file: {e}")
        return {}
    
    return charges


def get_unique_atom_types(molecules, selected_molnames=None):
    """Get unique atom types from selected molecules."""
    unique_types = []
    
    for mol_name, mol_data in molecules.items():
        if selected_molnames is None or mol_name in selected_molnames:
            for atom_data in mol_data['atoms'].values():
                atom_type = atom_data['type']
                # Skip special atoms like NETF, TORQ, and special marked atoms
                if (atom_type not in ['NETF', 'TORQ', 'M', 'MMM'] and 
                    not atom_data.get('special', False) and 
                    atom_type not in unique_types):
                    unique_types.append(atom_type)
    
    return unique_types


def get_element_from_atom_type(atom_type):
    """Extract element from atom type (first letter)."""
    return atom_type[0]


def generate_atomtypes_xml(atom_types):
    """Generate AtomTypes XML section."""
    xml_lines = ['<AtomTypes>']
    
    for atom_type in atom_types:
        element = get_element_from_atom_type(atom_type)
        xml_lines.append(f'<Type name="{atom_type}" class="{atom_type}" element="{element}" mass="0.0"/>')
    
    xml_lines.append('</AtomTypes>')
    return '\n'.join(xml_lines)


def create_unique_atom_names(molecules, selected_molnames=None):
    """Create unique atom names based on element with sequential numbering."""
    atom_id_to_name = {}
    element_counters = {}
    
    for mol_name, mol_data in molecules.items():
        if selected_molnames is None or mol_name in selected_molnames:
            # Sort atom IDs to ensure consistent ordering
            for atom_id in sorted(mol_data['atoms'].keys()):
                atom_data = mol_data['atoms'][atom_id]
                atom_type = atom_data['type']
                
                # Skip special atoms like NETF, TORQ, M, MMM, and special marked atoms
                if (atom_type in ['NETF', 'TORQ', 'M', 'MMM'] or 
                    atom_data.get('special', False)):
                    continue
                    
                element = atom_type[0]  # First letter is the element
                
                # Initialize counter for this element if not seen before
                if element not in element_counters:
                    element_counters[element] = 0
                
                # Create unique name
                unique_name = f"{element}{element_counters[element]}"
                
                # Store mapping
                atom_id_to_name[(mol_name, atom_id)] = unique_name
                
                # Increment counter for this element
                element_counters[element] += 1
    
    return atom_id_to_name


def generate_residues_xml(molecules, selected_molnames=None, molname_translations=None):
    """Generate Residues XML section without charges, using unique atom names."""
    xml_lines = ['<Residues>']
    
    # Create unique atom names
    atom_id_to_name = create_unique_atom_names(molecules, selected_molnames)
    
    for mol_name, mol_data in molecules.items():
        if selected_molnames is None or mol_name in selected_molnames:
            # Use translation if provided, otherwise use original name
            residue_name = mol_name
            if molname_translations and mol_name in molname_translations:
                residue_name = molname_translations[mol_name]
            xml_lines.append(f'<Residue name="{residue_name}">')
            
            # Add atoms (without charges, skip special atoms)
            for atom_id in sorted(mol_data['atoms'].keys()):
                atom_data = mol_data['atoms'][atom_id]
                atom_type = atom_data['type']
                
                # Skip special atoms
                if (atom_type in ['NETF', 'TORQ', 'M', 'MMM'] or 
                    atom_data.get('special', False)):
                    continue
                
                unique_name = atom_id_to_name[(mol_name, atom_id)]
                xml_lines.append(f'<Atom name="{unique_name}" type="{atom_type}"/>')
            
            # Add bonds (using unique names)
            for bond in mol_data['bonds']:
                atom1_id, atom2_id = bond['atoms']
                if ((mol_name, atom1_id) in atom_id_to_name and 
                    (mol_name, atom2_id) in atom_id_to_name):
                    atom1_name = atom_id_to_name[(mol_name, atom1_id)]
                    atom2_name = atom_id_to_name[(mol_name, atom2_id)]
                    xml_lines.append(f'<Bond atomName1="{atom1_name}" atomName2="{atom2_name}"/>')
            
            xml_lines.append('</Residue>')
    
    xml_lines.append('</Residues>')
    return '\n'.join(xml_lines)


def generate_nonbonded_force_xml(atom_types, charges):
    """Generate NonbondedForce XML section for charges."""
    xml_lines = ['<NonbondedForce coulomb14scale="1.0" lj14scale="1.0">']
    
    for atom_type in atom_types:
        charge = charges.get(atom_type, 0.0)
        xml_lines.append(f'<Atom type="{atom_type}" charge="{charge}" sigma="0.0" epsilon="0.0"/>')
    
    xml_lines.append('</NonbondedForce>')
    return '\n'.join(xml_lines)


def extract_bond_types(molecules, selected_molnames=None):
    """Extract unique bond types from molecule definitions."""
    bond_types = {}
    
    for mol_name, mol_data in molecules.items():
        if selected_molnames is None or mol_name in selected_molnames:
            for bond in mol_data['bonds']:
                atom1_id, atom2_id = bond['atoms']
                if atom1_id in mol_data['atoms'] and atom2_id in mol_data['atoms']:
                    atom1_type = mol_data['atoms'][atom1_id]['type']
                    atom2_type = mol_data['atoms'][atom2_id]['type']
                    
                    # Skip special atoms
                    if (atom1_type in ['NETF', 'TORQ', 'M', 'MMM'] or 
                        atom2_type in ['NETF', 'TORQ', 'M', 'MMM'] or
                        mol_data['atoms'][atom1_id].get('special', False) or
                        mol_data['atoms'][atom2_id].get('special', False)):
                        continue
                    
                    # Create normalized bond type key
                    sorted_types = tuple(sorted([atom1_type, atom2_type]))
                    if sorted_types not in bond_types:
                        bond_types[sorted_types] = {
                            'type1': sorted_types[0],
                            'type2': sorted_types[1]
                        }
    
    return bond_types


def extract_angle_types(molecules, selected_molnames=None):
    """Extract unique angle types from molecule definitions."""
    angle_types = {}
    
    for mol_name, mol_data in molecules.items():
        if selected_molnames is None or mol_name in selected_molnames:
            for angle in mol_data['angles']:
                atom1_id, atom2_id, atom3_id = angle['atoms']
                if (atom1_id in mol_data['atoms'] and 
                    atom2_id in mol_data['atoms'] and 
                    atom3_id in mol_data['atoms']):
                    
                    atom1_type = mol_data['atoms'][atom1_id]['type']
                    atom2_type = mol_data['atoms'][atom2_id]['type']
                    atom3_type = mol_data['atoms'][atom3_id]['type']
                    
                    # Skip special atoms
                    if any(t in ['NETF', 'TORQ', 'M', 'MMM'] for t in [atom1_type, atom2_type, atom3_type]):
                        continue
                    if any(mol_data['atoms'][aid].get('special', False) for aid in [atom1_id, atom2_id, atom3_id]):
                        continue
                    
                    # Create angle type key (center atom is always middle)
                    angle_key = f"{atom1_type}-{atom2_type}-{atom3_type}"
                    reverse_key = f"{atom3_type}-{atom2_type}-{atom1_type}"
                    
                    if angle_key not in angle_types and reverse_key not in angle_types:
                        angle_types[angle_key] = {
                            'type1': atom1_type,
                            'type2': atom2_type,
                            'type3': atom3_type
                        }
    
    return angle_types


def extract_dihedral_types(molecules, selected_molnames=None):
    """Extract unique dihedral types from molecule definitions."""
    dihedral_types = {}
    
    for mol_name, mol_data in molecules.items():
        if selected_molnames is None or mol_name in selected_molnames:
            for dihedral in mol_data['dihedrals']:
                atom1_id, atom2_id, atom3_id, atom4_id = dihedral['atoms']
                if all(aid in mol_data['atoms'] for aid in [atom1_id, atom2_id, atom3_id, atom4_id]):
                    
                    atom1_type = mol_data['atoms'][atom1_id]['type']
                    atom2_type = mol_data['atoms'][atom2_id]['type']
                    atom3_type = mol_data['atoms'][atom3_id]['type']
                    atom4_type = mol_data['atoms'][atom4_id]['type']
                    
                    # Skip special atoms
                    if any(t in ['NETF', 'TORQ', 'M', 'MMM'] for t in [atom1_type, atom2_type, atom3_type, atom4_type]):
                        continue
                    if any(mol_data['atoms'][aid].get('special', False) for aid in [atom1_id, atom2_id, atom3_id, atom4_id]):
                        continue
                    
                    # Create dihedral type key
                    dihedral_key = f"{atom1_type}-{atom2_type}-{atom3_type}-{atom4_type}"
                    reverse_key = f"{atom4_type}-{atom3_type}-{atom2_type}-{atom1_type}"
                    
                    if dihedral_key not in dihedral_types and reverse_key not in dihedral_types:
                        dihedral_types[dihedral_key] = {
                            'type1': atom1_type,
                            'type2': atom2_type,
                            'type3': atom3_type,
                            'type4': atom4_type
                        }
    
    return dihedral_types


def generate_harmonic_bond_force_xml(bond_types, parameters):
    """Generate HarmonicBondForce XML section."""
    xml_lines = ['<HarmonicBondForce>']
    
    for bond_key, bond_info in bond_types.items():
        type1, type2 = bond_info['type1'], bond_info['type2']
        
        # Look for parameter in bonds parameters
        param_key = f"{type1}_{type2}"
        if param_key not in parameters.get('bonds', {}):
            param_key = f"{type2}_{type1}"
        
        if param_key in parameters.get('bonds', {}):
            param_data = parameters['bonds'][param_key]
            if param_data['type'] == 'HAR' and len(param_data['values']) >= 2:
                length = param_data['values'][0]
                k = param_data['values'][1]
                xml_lines.append(f'<Bond type1="{type1}" type2="{type2}" length="{length}" k="{k}"/>')
    
    xml_lines.append('</HarmonicBondForce>')
    return '\n'.join(xml_lines)


def generate_harmonic_angle_force_xml(angle_types, parameters):
    """Generate HarmonicAngleForce XML section."""
    xml_lines = ['<HarmonicAngleForce>']
    
    for angle_key, angle_info in angle_types.items():
        type1, type2, type3 = angle_info['type1'], angle_info['type2'], angle_info['type3']
        
        # Look for parameter in angles parameters
        param_key = f"{type1}_{type2}_{type3}"
        if param_key not in parameters.get('angles', {}):
            param_key = f"{type3}_{type2}_{type1}"
        
        if param_key in parameters.get('angles', {}):
            param_data = parameters['angles'][param_key]
            if param_data['type'] == 'HAR' and len(param_data['values']) >= 2:
                theta_deg = float(param_data['values'][0])
                k = param_data['values'][1]
                
                # Convert angle from degrees to radians
                theta_rad = theta_deg * math.pi / 180.0
                xml_lines.append(f'<Angle class1="{type1}" class2="{type2}" class3="{type3}" angle="{theta_rad}" k="{k}"/>')
    
    xml_lines.append('</HarmonicAngleForce>')
    return '\n'.join(xml_lines)


def generate_periodic_torsion_force_xml(dihedral_types, parameters):
    """Generate PeriodicTorsionForce XML section."""
    xml_lines = ['<PeriodicTorsionForce>']
    
    for dihedral_key, dihedral_info in dihedral_types.items():
        type1, type2, type3, type4 = (dihedral_info['type1'], dihedral_info['type2'], 
                                      dihedral_info['type3'], dihedral_info['type4'])
        
        # Look for parameter in dihedrals parameters
        param_key = f"{type1}_{type2}_{type3}_{type4}"
        if param_key not in parameters.get('dihedrals', {}):
            param_key = f"{type4}_{type3}_{type2}_{type1}"
        
        if param_key in parameters.get('dihedrals', {}):
            param_data = parameters['dihedrals'][param_key]
            if param_data['type'] == 'NCO' and len(param_data['values']) >= 3:
                phase_deg = float(param_data['values'][0])
                k = param_data['values'][1]
                periodicity = param_data['values'][2]
                
                # Convert phase from degrees to radians
                phase_rad = phase_deg * math.pi / 180.0
                xml_lines.append(f'<Proper class1="{type1}" class2="{type2}" class3="{type3}" class4="{type4}" periodicity1="{periodicity}" phase1="{phase_rad}" k1="{k}"/>')
    
    xml_lines.append('</PeriodicTorsionForce>')
    return '\n'.join(xml_lines)


def create_interaction_matrices(nonbonded_interactions, interaction_type, atom_types):
    """Create symmetric matrices for nonbonded interactions."""
    n_types = len(atom_types)
    
    if interaction_type == 'SRD':
        # Group SRD parameters by power
        srd_by_power = defaultdict(list)
        for interaction in nonbonded_interactions['SRD']:
            if len(interaction['values']) >= 3:
                power = float(interaction['values'][1])  # P2: power
                srd_by_power[power].append(interaction)
        
        matrices_by_power = {}
        for power, params in srd_by_power.items():
            # Initialize matrices
            disp_matrix = [[0.0 for _ in range(n_types)] for _ in range(n_types)]
            r0_matrix = [[0.0 for _ in range(n_types)] for _ in range(n_types)]
            
            for param in params:
                atom1, atom2 = param['atom1'], param['atom2']
                if atom1 in atom_types and atom2 in atom_types:
                    i, j = atom_types.index(atom1), atom_types.index(atom2)
                    
                    # Convert units: kcal to kJ, Angstrom to nm
                    p1_value = float(param['values'][0]) * 4.184 * (0.1**(-power))
                    p3_value = float(param['values'][2]) * 0.1
                    
                    disp_matrix[i][j] = disp_matrix[j][i] = p1_value
                    r0_matrix[i][j] = r0_matrix[j][i] = p3_value
            
            matrices_by_power[power] = {'disp_matrix': disp_matrix, 'r0_matrix': r0_matrix}
        
        return matrices_by_power
    
    elif interaction_type == 'EXP':
        # Initialize matrices
        a_matrix = [[0.0 for _ in range(n_types)] for _ in range(n_types)]
        alpha_matrix = [[0.0 for _ in range(n_types)] for _ in range(n_types)]
        
        for interaction in nonbonded_interactions['EXP']:
            if len(interaction['values']) >= 2:
                atom1, atom2 = interaction['atom1'], interaction['atom2']
                if atom1 in atom_types and atom2 in atom_types:
                    i, j = atom_types.index(atom1), atom_types.index(atom2)
                    
                    # Convert units: kcal/mol to kJ/mol, Angstrom^-1 to nm^-1
                    a_value = float(interaction['values'][0]) * 4.184
                    alpha_value = float(interaction['values'][1]) / 0.1
                    
                    a_matrix[i][j] = a_matrix[j][i] = a_value
                    alpha_matrix[i][j] = alpha_matrix[j][i] = alpha_value
        
        return {'a_matrix': a_matrix, 'alpha_matrix': alpha_matrix}


def matrix_to_xml_string(matrix):
    """Convert matrix to XML string format (flattened)."""
    result = []
    for row in matrix:
        result.extend([str(val) for val in row])
    return '\t'.join(result)


def generate_srd_custom_force_xml(power, matrices, atom_types):
    """Generate CustomNonbondedForce XML for SRD interactions."""
    n_types = len(atom_types)
    disp_matrix = matrices['disp_matrix']
    r0_matrix = matrices['r0_matrix']
    
    power_str = str(int(power)) if power == int(power) else str(power)
    abs_power_str = str(int(abs(power))) if abs(power) == int(abs(power)) else str(abs(power))
    
    xml_lines = []
    xml_lines.append(f'<CustomNonbondedForce energy="disp/(r^{abs_power_str} + r0^{abs_power_str}); disp=dispTable(t1,t2); r0=rTable(t1,t2)" bondCutoff="2">')
    
    # Dispersion table
    xml_lines.append(f'<Function name="dispTable" type="Discrete2D" xsize="{n_types}" ysize="{n_types}">')
    xml_lines.append(matrix_to_xml_string(disp_matrix))
    xml_lines.append('</Function>')
    
    # R0 table
    xml_lines.append(f'<Function name="rTable" type="Discrete2D" xsize="{n_types}" ysize="{n_types}">')
    xml_lines.append(matrix_to_xml_string(r0_matrix))
    xml_lines.append('</Function>')
    
    # Per-particle parameter
    xml_lines.append('<PerParticleParameter name="t"/>')
    
    # Atom type assignments
    for i, atom_type in enumerate(atom_types):
        xml_lines.append(f'<Atom type="{atom_type}" t="{i}"/>')
    
    xml_lines.append('</CustomNonbondedForce>')
    return '\n'.join(xml_lines)


def generate_exp_custom_force_xml(matrices, atom_types):
    """Generate CustomNonbondedForce XML for EXP interactions."""
    n_types = len(atom_types)
    a_matrix = matrices['a_matrix']
    alpha_matrix = matrices['alpha_matrix']
    
    xml_lines = []
    xml_lines.append('<CustomNonbondedForce energy="a*exp(-alpha*r); a=aTable(t1,t2); alpha=alphaTable(t1,t2)" bondCutoff="2">')
    
    # A table (amplitude)
    xml_lines.append(f'<Function name="aTable" type="Discrete2D" xsize="{n_types}" ysize="{n_types}">')
    xml_lines.append(matrix_to_xml_string(a_matrix))
    xml_lines.append('</Function>')
    
    # Alpha table (decay parameter)
    xml_lines.append(f'<Function name="alphaTable" type="Discrete2D" xsize="{n_types}" ysize="{n_types}">')
    xml_lines.append(matrix_to_xml_string(alpha_matrix))
    xml_lines.append('</Function>')
    
    # Per-particle parameter
    xml_lines.append('<PerParticleParameter name="t"/>')
    
    # Atom type assignments
    for i, atom_type in enumerate(atom_types):
        xml_lines.append(f'<Atom type="{atom_type}" t="{i}"/>')
    
    xml_lines.append('</CustomNonbondedForce>')
    return '\n'.join(xml_lines)


def main():
    parser = argparse.ArgumentParser(
        description="Convert CRYOFF .off force field files to OpenMM .xml format"
    )
    parser.add_argument("-off", required=True, help="Input .off file path")
    parser.add_argument("-output", required=True, help="Output .xml file path")
    parser.add_argument(
        "-molnames", 
        help="Comma-separated list of molecule names to include (e.g., 'UNK,CYCQM')"
    )
    parser.add_argument(
        "-charges", 
        help="File with Atom Name and Charge columns for charge assignment"
    )
    parser.add_argument(
        "-molname_translations",
        help="Comma-separated list of 3-letter residue names corresponding to molnames (e.g., 'CYC,SOL')"
    )
    
    args = parser.parse_args()
    
    # Parse selected molecule names
    selected_molnames = None
    molname_translations = {}
    if args.molnames:
        selected_molnames = [name.strip() for name in args.molnames.split(',')]
        
        # Handle molname translations
        if args.molname_translations:
            translations = [name.strip() for name in args.molname_translations.split(',')]
            if len(translations) != len(selected_molnames):
                print(f"Error: Number of translations ({len(translations)}) must match number of molnames ({len(selected_molnames)})")
                sys.exit(1)
            molname_translations = dict(zip(selected_molnames, translations))
        else:
            # Use first 3 letters of each molname as default translation
            molname_translations = {name: name[:3] for name in selected_molnames}
    
    # Parse .off file
    print(f"Parsing {args.off}...")
    parser = OffFileParser(args.off)
    molecules, parameters, nonbonded_interactions = parser.parse()
    
    print(f"Found {len(molecules)} molecule types: {list(molecules.keys())}")
    
    # Read charges if provided
    charges = {}
    if args.charges:
        print(f"Reading charges from {args.charges}...")
        charges = read_charges_file(args.charges)
        print(f"Loaded charges for {len(charges)} atom types")
    
    # Get unique atom types from selected molecules
    atom_types = get_unique_atom_types(molecules, selected_molnames)
    print(f"Found {len(atom_types)} unique atom types: {atom_types}")
    
    # Generate XML sections
    xml_sections = []
    
    # AtomTypes section
    xml_sections.append(generate_atomtypes_xml(atom_types))
    
    # Residues section (without charges)
    xml_sections.append(generate_residues_xml(molecules, selected_molnames, molname_translations))
    
    # NonbondedForce section for charges
    if charges:
        xml_sections.append(generate_nonbonded_force_xml(atom_types, charges))
    
    # Extract bonded interaction types
    bond_types = extract_bond_types(molecules, selected_molnames)
    angle_types = extract_angle_types(molecules, selected_molnames)
    dihedral_types = extract_dihedral_types(molecules, selected_molnames)
    
    # Generate bonded force sections
    if bond_types and parameters.get('bonds'):
        xml_sections.append(generate_harmonic_bond_force_xml(bond_types, parameters))
    
    if angle_types and parameters.get('angles'):
        xml_sections.append(generate_harmonic_angle_force_xml(angle_types, parameters))
    
    if dihedral_types and parameters.get('dihedrals'):
        xml_sections.append(generate_periodic_torsion_force_xml(dihedral_types, parameters))
    
    # Generate CustomNonbondedForce sections for EXP and SRD
    if nonbonded_interactions['EXP']:
        print("Generating EXP CustomNonbondedForce...")
        exp_matrices = create_interaction_matrices(nonbonded_interactions, 'EXP', atom_types)
        xml_sections.append(generate_exp_custom_force_xml(exp_matrices, atom_types))
    
    if nonbonded_interactions['SRD']:
        print("Generating SRD CustomNonbondedForce...")
        srd_matrices_by_power = create_interaction_matrices(nonbonded_interactions, 'SRD', atom_types)
        for power in sorted(srd_matrices_by_power.keys()):
            xml_sections.append(generate_srd_custom_force_xml(power, srd_matrices_by_power[power], atom_types))
    
    # Write output file
    print(f"Writing {args.output}...")
    with open(args.output, 'w') as f:
        f.write('<ForceField>\n\n')
        f.write('\n\n'.join(xml_sections) + '\n')
        f.write('\n</ForceField>\n')
    
    print("Conversion complete!")


if __name__ == "__main__":
    main()