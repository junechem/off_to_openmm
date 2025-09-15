#!/usr/bin/env python3
"""
Convert CRYOFF .off force field files to OpenMM .xml format.
"""

import argparse
import sys
import re
import math
from collections import defaultdict

# Try to import OpenMM for atomic masses
try:
    from openmm.app import Element
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False
    print("Warning: OpenMM not found. Atom masses will be set to 0.0.")


class OffFileParser:
    """Parser for CRYOFF .off files with support for multiple format variations."""
    
    def __init__(self, filename):
        self.filename = filename
        self.molecules = {}
        self.parameters = {}
        self.nonbonded_interactions = {'COU': [], 'EXP': [], 'SRD': [], 'BUC': [], 'STR': [], 'POW': []}
        
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
            if re.search(r'\[\s*ATO', line):
                i = self._parse_atoms_section(lines, i + 1, mol_name)
            # Parse bonds section
            elif re.search(r'\[\s*BON', line):
                i = self._parse_bonds_section(lines, i + 1, mol_name)
            # Parse angles section
            elif re.search(r'\[\s*ANG', line):
                i = self._parse_angles_section(lines, i + 1, mol_name)
            # Parse dihedrals section
            elif re.search(r'\[\s*DIH', line):
                i = self._parse_dihedrals_section(lines, i + 1, mol_name)
            # Parse exclusions section
            elif re.search(r'\[\s*EXC', line):
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
                    is_virtual_site = '*' in parts[0]
                    
                    atom_data = {
                        'name': atom_name,
                        'type': atom_type,
                        'special': is_virtual_site
                    }
                    
                    # Parse virtual site definition if present
                    if is_virtual_site and len(parts) > 3:
                        vsite_def = self._parse_virtual_site_definition(' '.join(parts[3:]))
                        if vsite_def:
                            atom_data['virtual_site'] = vsite_def
                    
                    self.molecules[mol_name]['atoms'][atom_id] = atom_data
                except ValueError:
                    # Skip lines that can't be parsed as atoms
                    pass
            i += 1
            
        return i
    
    def _parse_virtual_site_definition(self, definition):
        """Parse virtual site definition from .off format.
        
        Format: "3: 0.6 1 0.2 2 0.2 3" means 3 atoms with weights 0.6, 0.2, 0.2
        for atoms 1, 2, 3 respectively.
        """
        try:
            parts = definition.split()
            if len(parts) < 2:
                return None
                
            # Extract number of atoms and colon
            num_atoms_str = parts[0].rstrip(':')
            num_atoms = int(num_atoms_str)
            
            # Parse weight-atom pairs
            weights = []
            atom_refs = []
            
            i = 1
            while i < len(parts) and len(weights) < num_atoms:
                if i < len(parts):
                    weight = float(parts[i])
                    weights.append(weight)
                if i + 1 < len(parts):
                    atom_ref = int(parts[i + 1])
                    atom_refs.append(atom_ref)
                i += 2
            
            if len(weights) == num_atoms and len(atom_refs) == num_atoms:
                return {
                    'num_atoms': num_atoms,
                    'weights': weights,
                    'atom_refs': atom_refs
                }
                
        except (ValueError, IndexError):
            pass
            
        return None
    
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
                
            if line.startswith(('HAR', 'QUA')):
                # New bond group definition (HAR = harmonic, QUA = quartic)
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
        """Parse nonbonded interactions from Inter-Potential section at bottom of file."""
        lines = content.split('\n')
        
        # Find the Inter-Potential section
        start_idx = -1
        for i, line in enumerate(lines):
            if line.strip().startswith('Inter-Potential:'):
                start_idx = i + 1
                break
        
        if start_idx == -1:
            print("Warning: Inter-Potential section not found in .off file")
            return
        
        # Parse interactions until we hit Molecular-Definition or end of file
        for i in range(start_idx, len(lines)):
            line = lines[i].strip()
            
            # Stop at Molecular-Definition or empty sections
            if line.startswith('Molecular-Definition:') or not line:
                continue
            
            # Parse interaction lines: ATOM1~ATOM2: TYPE param1 param2 param3 Min: X Max: Y
            if '~' in line and ':' in line:
                self._parse_inter_potential_line(line)
    
    def _parse_inter_potential_line(self, line):
        """Parse a single line from Inter-Potential section.
        
        Format: ATOM1~ATOM2: TYPE param1 param2 param3 Min: X Max: Y
        """
        try:
            # Split at the colon to separate atom pair from the rest
            atom_part, rest = line.split(':', 1)
            atom_pair = atom_part.strip()
            
            # Extract atom types from ATOM1~ATOM2 format
            if '~' not in atom_pair:
                return
            atom1, atom2 = atom_pair.split('~', 1)
            
            # Parse the rest: TYPE param1 param2 param3 Min: X Max: Y
            parts = rest.strip().split()
            if len(parts) < 2:  # Need at least interaction type and one parameter
                return
                
            interaction_type = parts[0]
            
            # Find where "Min:" appears to know where parameters end
            min_idx = -1
            for i, part in enumerate(parts):
                if part == 'Min:':
                    min_idx = i
                    break
            
            # Extract parameters (everything between interaction type and Min:)
            if min_idx > 1:
                values = parts[1:min_idx]
            else:
                # If no Min: found, take all remaining parts as parameters
                values = parts[1:]
            
            # Map interaction type to our categories
            if interaction_type == 'COU':
                interaction_key = 'COU'
            elif interaction_type == 'EXP':
                interaction_key = 'EXP'
            elif interaction_type == 'SRD':
                interaction_key = 'SRD'
            elif interaction_type in ['BUC', 'BUCK']:
                interaction_key = 'BUC'
            elif interaction_type in ['STR', 'STRC']:
                interaction_key = 'STR'
            elif interaction_type == 'POW':
                interaction_key = 'POW'
            else:
                # Skip unknown interaction types
                return
            
            # Add to nonbonded interactions
            self.nonbonded_interactions[interaction_key].append({
                'atom1': atom1,
                'atom2': atom2,
                'fix_flag': 'FIX',  # Not specified in Inter-Potential format
                'values': values
            })
            
        except (ValueError, IndexError) as e:
            # Skip malformed lines
            pass


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
                # Skip only NETF and TORQ atoms, but include virtual sites (M, MMM)
                if (atom_type not in ['NETF', 'TORQ'] and 
                    atom_type not in unique_types):
                    unique_types.append(atom_type)
    
    return unique_types


def get_element_from_atom_type(atom_type):
    """Extract element from atom type (first letter)."""
    return atom_type[0]


def get_atomic_mass(element_symbol):
    """Get atomic mass for an element symbol using OpenMM's Element class."""
    if not OPENMM_AVAILABLE:
        return 0.0
    
    try:
        element = Element.getBySymbol(element_symbol)
        # element.mass returns a Quantity object, need to get the value in atomic mass units
        return float(element.mass.value_in_unit(element.mass.unit))
    except Exception:
        # Element not found or other error, return 0.0 as fallback
        return 0.0


def is_virtual_site_type(atom_type, molecules, selected_molnames=None):
    """Check if an atom type belongs to virtual sites."""
    for mol_name, mol_data in molecules.items():
        if selected_molnames is None or mol_name in selected_molnames:
            for atom_data in mol_data['atoms'].values():
                if atom_data['type'] == atom_type:
                    # Check if this atom is a virtual site
                    if atom_data.get('special') or atom_data.get('virtual_site'):
                        return True
    return False


def generate_atomtypes_xml(atom_types, molecules, selected_molnames=None):
    """Generate AtomTypes XML section."""
    xml_lines = ['<AtomTypes>']
    
    for atom_type in atom_types:
        # Virtual sites should not have an element attribute and should have mass="0.0"
        if is_virtual_site_type(atom_type, molecules, selected_molnames):
            xml_lines.append(f'<Type name="{atom_type}" class="{atom_type}" mass="0.0"/>')
        else:
            element = get_element_from_atom_type(atom_type)
            mass = get_atomic_mass(element)
            xml_lines.append(f'<Type name="{atom_type}" class="{atom_type}" element="{element}" mass="{mass}"/>')
    
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
                
                # Skip only NETF and TORQ atoms, but include virtual sites (M, MMM)
                if atom_type in ['NETF', 'TORQ']:
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
            
            # Add atoms (without charges, skip only NETF/TORQ atoms)
            virtual_sites = []  # Collect virtual sites for later processing
            
            for atom_id in sorted(mol_data['atoms'].keys()):
                atom_data = mol_data['atoms'][atom_id]
                atom_type = atom_data['type']
                
                # Skip only NETF and TORQ atoms, but include virtual sites (M, MMM)
                if atom_type in ['NETF', 'TORQ']:
                    continue
                
                unique_name = atom_id_to_name[(mol_name, atom_id)]
                xml_lines.append(f'<Atom name="{unique_name}" type="{atom_type}"/>')
                
                # Collect virtual site information
                if atom_data.get('virtual_site'):
                    virtual_sites.append((atom_id, atom_data, unique_name))
            
            # Add bonds (using unique names)
            for bond in mol_data['bonds']:
                atom1_id, atom2_id = bond['atoms']
                if ((mol_name, atom1_id) in atom_id_to_name and 
                    (mol_name, atom2_id) in atom_id_to_name):
                    atom1_name = atom_id_to_name[(mol_name, atom1_id)]
                    atom2_name = atom_id_to_name[(mol_name, atom2_id)]
                    xml_lines.append(f'<Bond atomName1="{atom1_name}" atomName2="{atom2_name}"/>')
            
            # Add virtual sites
            for vsite_atom_id, vsite_data, vsite_name in virtual_sites:
                vsite_def = vsite_data['virtual_site']
                
                # Convert atom references to atom names
                ref_atom_names = []
                valid_vsite = True
                
                for ref_atom_id in vsite_def['atom_refs']:
                    if (mol_name, ref_atom_id) in atom_id_to_name:
                        ref_atom_names.append(atom_id_to_name[(mol_name, ref_atom_id)])
                    else:
                        valid_vsite = False
                        break
                
                if valid_vsite and vsite_def['num_atoms'] == 3:
                    # Generate average3 virtual site
                    xml_lines.append(
                        f'<VirtualSite type="average3" siteName="{vsite_name}" '
                        f'atomName1="{ref_atom_names[0]}" '
                        f'atomName2="{ref_atom_names[1]}" '
                        f'atomName3="{ref_atom_names[2]}" '
                        f'weight1="{vsite_def["weights"][0]}" '
                        f'weight2="{vsite_def["weights"][1]}" '
                        f'weight3="{vsite_def["weights"][2]}"/>'
                    )
                elif valid_vsite and vsite_def['num_atoms'] == 2:
                    # Generate average2 virtual site
                    xml_lines.append(
                        f'<VirtualSite type="average2" siteName="{vsite_name}" '
                        f'atomName1="{ref_atom_names[0]}" '
                        f'atomName2="{ref_atom_names[1]}" '
                        f'weight1="{vsite_def["weights"][0]}" '
                        f'weight2="{vsite_def["weights"][1]}"/>'
                    )
            
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


def generate_custom_bond_force_xml(bond_types, parameters):
    """Generate CustomBondForce XML section for QUA bonds."""
    qua_bonds = []
    
    # First, check if there are any QUA bonds
    for bond_key, bond_info in bond_types.items():
        type1, type2 = bond_info['type1'], bond_info['type2']
        
        # Look for QUA parameter in bonds parameters
        param_key = f"{type1}_{type2}"
        if param_key not in parameters.get('bonds', {}):
            param_key = f"{type2}_{type1}"
        
        if param_key in parameters.get('bonds', {}):
            param_data = parameters['bonds'][param_key]
            if param_data['type'] == 'QUA' and len(param_data['values']) >= 4:
                # QUA parameters from #define section (already in correct units)
                r0 = float(param_data['values'][0])  # nm
                k2 = float(param_data['values'][1])  # kJ/(mol*nm²)
                k3 = float(param_data['values'][2])  # kJ/(mol*nm³)
                k4 = float(param_data['values'][3])  # kJ/(mol*nm⁴)
                
                qua_bonds.append(f'<Bond type1="{type1}" type2="{type2}" r0="{r0}" k2="{k2}" k3="{k3}" k4="{k4}"/>')
    
    # Only generate CustomBondForce if there are QUA bonds
    if not qua_bonds:
        return ""
    
    xml_lines = ['<CustomBondForce energy="(k2/2)*(r-r0)^2 + (k3/3)*(r-r0)^3 + (k4/4)*(r-r0)^4">']
    xml_lines.append('<PerBondParameter name="r0"/>')
    xml_lines.append('<PerBondParameter name="k2"/>')
    xml_lines.append('<PerBondParameter name="k3"/>')
    xml_lines.append('<PerBondParameter name="k4"/>')
    xml_lines.extend(qua_bonds)
    xml_lines.append('</CustomBondForce>')
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


def process_buc_interactions(nonbonded_interactions):
    """Process BUC interactions and add equivalent EXP and SRD entries.
    
    BUC interaction: U(r) = P2/r^6 + P1*exp(-P3*r)
    - P1*exp(-P3*r) goes to EXP section
    - P2/r^6 goes to SRD section with r0=0
    """
    for interaction in nonbonded_interactions['BUC']:
        if len(interaction['values']) >= 3:
            atom1 = interaction['atom1']
            atom2 = interaction['atom2']
            fix_flag = interaction['fix_flag']
            
            # BUC parameters: P1, P2, P3
            p1 = float(interaction['values'][0])  # kcal/mol (for exp term)
            p2 = float(interaction['values'][1])  # kcal*angstrom^6/mol (for r^-6 term)
            p3 = float(interaction['values'][2])  # 1/angstrom (for exp term)
            
            # Add P1*exp(-P3*r) to EXP section
            exp_entry = {
                'atom1': atom1,
                'atom2': atom2,
                'fix_flag': fix_flag,
                'values': [str(p1), str(p3)]
            }
            nonbonded_interactions['EXP'].append(exp_entry)
            
            # Add P2/r^6 to SRD section with r0=0
            srd_entry = {
                'atom1': atom1,
                'atom2': atom2,
                'fix_flag': fix_flag,
                'values': [str(p2), '-6', '0.0']  # P1, power, r0
            }
            nonbonded_interactions['SRD'].append(srd_entry)


def process_pow_interactions(nonbonded_interactions):
    """Process POW interactions and add equivalent SRD entries.
    
    POW interaction: U(r) = P1*r^P2
    - This goes to SRD section with r0=0 as P1/(r^-P2 + 0^-P2) = P1*r^P2
    """
    for interaction in nonbonded_interactions['POW']:
        if len(interaction['values']) >= 2:
            atom1 = interaction['atom1']
            atom2 = interaction['atom2']
            fix_flag = interaction['fix_flag']
            
            # POW parameters: P1, P2
            p1 = float(interaction['values'][0])  # kcal*angstrom^P2/mol 
            p2 = float(interaction['values'][1])  # dimensionless power
            
            # POW: U(r) = P1*r^P2
            # Convert to SRD format: P1/(r^power + r0^power) with r0=0
            # SRD unit conversion will handle: P1 * 4.184 * (0.1**(-power))
            # We want: P1 * 4.184 * (0.1**abs(P2))
            # So we need: power = -abs(P2) so that (0.1**(-power)) = (0.1**abs(P2))
            srd_entry = {
                'atom1': atom1,
                'atom2': atom2,
                'fix_flag': fix_flag,
                'values': [str(p1), str(p2), '0.0']  # P1, P2 (as-is), r0=0
            }
            nonbonded_interactions['SRD'].append(srd_entry)


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
    
    elif interaction_type == 'STR':
        # Initialize matrices
        p1_matrix = [[0.0 for _ in range(n_types)] for _ in range(n_types)]
        p2_matrix = [[0.0 for _ in range(n_types)] for _ in range(n_types)]
        p3_matrix = [[0.0 for _ in range(n_types)] for _ in range(n_types)]
        
        for interaction in nonbonded_interactions['STR']:
            if len(interaction['values']) >= 3:
                atom1, atom2 = interaction['atom1'], interaction['atom2']
                if atom1 in atom_types and atom2 in atom_types:
                    i, j = atom_types.index(atom1), atom_types.index(atom2)
                    
                    # Unit conversions for STR parameters
                    p1_off = float(interaction['values'][0])  # kcal/mol (in Angstrom^P2 units)
                    p2_value = float(interaction['values'][1])  # dimensionless
                    p3_off = float(interaction['values'][2])   # Angstrom
                    
                    # Convert P1: kcal/mol*Angstrom^P2 to kJ/mol*nm^P2
                    p1_value = p1_off * 4.184 * (0.1**p2_value)
                    # P2 remains dimensionless
                    # Convert P3: Angstrom to nm
                    p3_value = p3_off * 0.1
                    
                    p1_matrix[i][j] = p1_matrix[j][i] = p1_value
                    p2_matrix[i][j] = p2_matrix[j][i] = p2_value
                    p3_matrix[i][j] = p3_matrix[j][i] = p3_value
        
        return {'p1_matrix': p1_matrix, 'p2_matrix': p2_matrix, 'p3_matrix': p3_matrix}


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


def generate_str_custom_force_xml(matrices, atom_types):
    """Generate CustomNonbondedForce XML for STR interactions.
    
    STR interaction:
    if r <= P3: U(r) = P1*((1/r^P2) - (1/P3^P2) + (P2*(r - P3))/(P3^(P2 + 1)))
    if r > P3:  U(r) = 0
    """
    n_types = len(atom_types)
    p1_matrix = matrices['p1_matrix']
    p2_matrix = matrices['p2_matrix'] 
    p3_matrix = matrices['p3_matrix']
    
    xml_lines = []
    # Use step function: step(P3 - r) = 1 if r <= P3, 0 if r > P3
    energy_expr = "step(p3 - r) * p1 * ((1/r^p2) - (1/p3^p2) + (p2*(r - p3))/(p3^(p2 + 1))); "
    energy_expr += "p1=p1Table(t1,t2); p2=p2Table(t1,t2); p3=p3Table(t1,t2)"
    
    xml_lines.append(f'<CustomNonbondedForce energy="{energy_expr}" bondCutoff="2">')
    
    # P1 table
    xml_lines.append(f'<Function name="p1Table" type="Discrete2D" xsize="{n_types}" ysize="{n_types}">')
    xml_lines.append(matrix_to_xml_string(p1_matrix))
    xml_lines.append('</Function>')
    
    # P2 table  
    xml_lines.append(f'<Function name="p2Table" type="Discrete2D" xsize="{n_types}" ysize="{n_types}">')
    xml_lines.append(matrix_to_xml_string(p2_matrix))
    xml_lines.append('</Function>')
    
    # P3 table
    xml_lines.append(f'<Function name="p3Table" type="Discrete2D" xsize="{n_types}" ysize="{n_types}">')
    xml_lines.append(matrix_to_xml_string(p3_matrix))
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
    xml_sections.append(generate_atomtypes_xml(atom_types, molecules, selected_molnames))
    
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
        custom_bond_xml = generate_custom_bond_force_xml(bond_types, parameters)
        if custom_bond_xml:  # Only add if there are QUA bonds
            xml_sections.append(custom_bond_xml)
    
    if angle_types and parameters.get('angles'):
        xml_sections.append(generate_harmonic_angle_force_xml(angle_types, parameters))
    
    if dihedral_types and parameters.get('dihedrals'):
        xml_sections.append(generate_periodic_torsion_force_xml(dihedral_types, parameters))
    
    # Process BUC interactions and convert to EXP and SRD components
    if nonbonded_interactions['BUC']:
        print("Processing BUC interactions...")
        process_buc_interactions(nonbonded_interactions)
    
    # Process POW interactions and convert to SRD components  
    if nonbonded_interactions['POW']:
        print("Processing POW interactions...")
        process_pow_interactions(nonbonded_interactions)
    
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
    
    if nonbonded_interactions['STR']:
        print("Generating STR CustomNonbondedForce...")
        str_matrices = create_interaction_matrices(nonbonded_interactions, 'STR', atom_types)
        xml_sections.append(generate_str_custom_force_xml(str_matrices, atom_types))
    
    # Write output file
    print(f"Writing {args.output}...")
    with open(args.output, 'w') as f:
        f.write('<ForceField>\n\n')
        f.write('\n\n'.join(xml_sections) + '\n')
        f.write('\n</ForceField>\n')
    
    print("Conversion complete!")


if __name__ == "__main__":
    main()