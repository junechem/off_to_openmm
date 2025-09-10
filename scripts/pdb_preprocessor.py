#!/usr/bin/env python3
"""
PDB Preprocessor for OpenMM Compatibility

This script processes PDB files to:
1. Convert atom classes to atom types if needed
2. Add CONECT records based on force field XML definitions

Usage:
    python pdb_preprocessor.py -pdb input.pdb -xml forcefield.xml -output output.pdb
"""

import argparse
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


class PDBPreprocessor:
    """Preprocesses PDB files for OpenMM compatibility."""
    
    def __init__(self, xml_file):
        """Initialize with force field XML file."""
        self.xml_file = xml_file
        self.atom_types = {}  # class -> type mapping
        self.residue_bonds = {}  # residue_name -> [(atom1, atom2), ...]
        self.residue_atom_types = {}  # residue_name -> {atom_name: atom_type}
        self.type_to_atoms = {}  # residue_name -> {atom_type: [atom_names]}
        self._parse_xml()
    
    def _parse_xml(self):
        """Parse the force field XML file to extract atom types and bonds."""
        try:
            tree = ET.parse(self.xml_file)
            root = tree.getroot()
            
            # Parse AtomTypes to get class -> type mapping
            for atom_type in root.findall('.//AtomTypes/Type'):
                class_name = atom_type.get('class')
                type_name = atom_type.get('name')
                if class_name and type_name:
                    self.atom_types[class_name] = type_name
            
            # Parse Residues to get bonds and atom type mapping for each residue
            for residue in root.findall('.//Residues/Residue'):
                residue_name = residue.get('name')
                if not residue_name:
                    continue
                
                # Get all atoms in this residue and their types
                atom_types_map = {}
                type_to_atoms_map = {}
                
                for atom in residue.findall('Atom'):
                    atom_name = atom.get('name')
                    atom_type = atom.get('type')
                    if atom_name and atom_type:
                        atom_types_map[atom_name] = atom_type
                        
                        # Build reverse mapping: type -> [atom_names]
                        if atom_type not in type_to_atoms_map:
                            type_to_atoms_map[atom_type] = []
                        type_to_atoms_map[atom_type].append(atom_name)
                
                # Get all bonds in this residue
                bonds = []
                for bond in residue.findall('Bond'):
                    atom1 = bond.get('atomName1')
                    atom2 = bond.get('atomName2')
                    if atom1 and atom2:
                        bonds.append((atom1, atom2))
                
                self.residue_atom_types[residue_name] = atom_types_map
                self.type_to_atoms[residue_name] = type_to_atoms_map
                self.residue_bonds[residue_name] = bonds
                
        except ET.ParseError as e:
            print(f"Error parsing XML file: {e}")
            sys.exit(1)
        except FileNotFoundError:
            print(f"XML file not found: {self.xml_file}")
            sys.exit(1)
    
    def _is_atom_class(self, atom_name):
        """Check if an atom name appears to be a class rather than type."""
        # If the atom name matches a known class, it's likely a class
        return atom_name in self.atom_types
    
    def _convert_atom_name(self, atom_name):
        """Convert atom class to atom type if needed."""
        if atom_name in self.atom_types:
            return self.atom_types[atom_name]
        return atom_name
    
    def _find_xml_atom_name_for_type(self, residue_name, atom_type):
        """Find XML atom name for a given atom type in a residue."""
        if residue_name not in self.type_to_atoms:
            return None
        
        if atom_type in self.type_to_atoms[residue_name]:
            # Return the first XML atom name of this type
            xml_atom_names = self.type_to_atoms[residue_name][atom_type]
            if xml_atom_names:
                return xml_atom_names[0]
        
        return None
    
    def _find_xml_atom_name_for_type_with_counter(self, residue_name, atom_type, type_counters):
        """Find XML atom name for a given atom type, using counters to track multiple atoms of same type."""
        if residue_name not in self.type_to_atoms:
            return None
        
        if atom_type in self.type_to_atoms[residue_name]:
            xml_atom_names = self.type_to_atoms[residue_name][atom_type]
            if xml_atom_names:
                # Get current count for this type
                current_count = type_counters.get(atom_type, 0)
                
                # If we have enough XML atom names, return the next one
                if current_count < len(xml_atom_names):
                    xml_atom_name = xml_atom_names[current_count]
                    type_counters[atom_type] = current_count + 1
                    return xml_atom_name
        
        return None
    
    def process_pdb(self, input_pdb, output_pdb):
        """Process PDB file: convert classes to types and add CONECT records."""
        atom_records = []
        other_records = []
        atom_index_map = {}  # (residue, atom_name) -> pdb_index
        
        # Read and parse PDB file
        try:
            with open(input_pdb, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            print(f"PDB file not found: {input_pdb}")
            sys.exit(1)
        
        # Process each line with simple approach
        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                # Parse ATOM/HETATM record (PDB format specification)
                record_type = line[0:6].strip()
                atom_index = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()
                chain_id = line[21:22]
                residue_seq = line[22:26].strip()
                
                original_atom_name = atom_name
                converted_atom_name = atom_name
                
                # Check if atom name is a class and convert if needed
                if self._is_atom_class(atom_name):
                    converted_atom_name = self._convert_atom_name(atom_name)
                    # Update the line with new atom name
                    line = line[:12] + f"{converted_atom_name:>4}" + line[16:]
                
                atom_records.append({
                    'line': line,
                    'index': atom_index,
                    'name': converted_atom_name,
                    'original_name': original_atom_name,
                    'residue': residue_name,
                    'chain': chain_id,
                    'seq': residue_seq
                })
                
                # Map for CONECT generation
                key = (residue_name, atom_name)
                atom_index_map[key] = atom_index
                
            elif not line.startswith(('TER', 'CONECT', 'END')):
                # Keep other records but skip TER, CONECT, END since we'll regenerate them
                other_records.append(line)
        
        
        # Generate CONECT records
        conect_records = self._generate_conect_records(atom_records, atom_index_map)
        
        # Write output PDB
        try:
            with open(output_pdb, 'w') as f:
                # Write header records first (TITLE, REMARK, CRYST1, MODEL, etc.)
                for record in other_records:
                    if not record.startswith('END'):
                        f.write(record)
                
                # Write ATOM/HETATM records
                for atom in atom_records:
                    f.write(atom['line'])
                
                # Write CONECT records
                for conect in conect_records:
                    f.write(conect)
                
                # Add END record
                f.write("END\n")
                
        except IOError as e:
            print(f"Error writing output file: {e}")
            sys.exit(1)
        
        print(f"Processed PDB file: {input_pdb} -> {output_pdb}")
        print(f"Added {len(conect_records)} CONECT records")
        
        # Report any atom class conversions
        converted_atoms = [atom for atom in atom_records 
                          if atom['name'] != atom['original_name']]
        if converted_atoms:
            print(f"Converted {len(converted_atoms)} atom classes to types:")
            for atom in converted_atoms[:5]:  # Show first 5
                print(f"  {atom['original_name']} -> {atom['name']}")
            if len(converted_atoms) > 5:
                print(f"  ... and {len(converted_atoms) - 5} more")
    
    def _generate_conect_records(self, atom_records, atom_index_map):
        """Generate CONECT records based on residue bond definitions."""
        conect_records = []
        bonds_added = set()  # Track bonds globally to avoid duplicates
        
        # Get unique residue types present in PDB
        unique_residues = set(atom['residue'] for atom in atom_records)
        xml_residues = set(self.residue_bonds.keys())
        
        # Check for residue name mismatches and provide warnings
        missing_residues = unique_residues - xml_residues
        if missing_residues:
            print(f"Warning: PDB contains residues not found in XML: {', '.join(missing_residues)}")
            print("Available XML residues:", ', '.join(xml_residues))
            print("Continuing with available residue definitions...")
        
        # Group atoms by residue instance
        residues = {}
        for atom in atom_records:
            res_key = (atom['residue'], atom['chain'], atom['seq'])
            if res_key not in residues:
                residues[res_key] = []
            residues[res_key].append(atom)
        
        # Generate CONECT records for each residue instance
        total_bonds_processed = 0
        for (residue_name, chain, seq), atoms in residues.items():
            if residue_name not in self.residue_bonds:
                continue  # Skip residues we don't have definitions for
            
            # Create mapping from XML atom names to PDB indices for this residue instance
            xml_atom_to_pdb_index = {}
            
            # Now that atom names have been converted to XML format during processing,
            # we can directly map XML atom names to PDB indices
            for atom in atoms:
                atom_name = atom['name']  # This is now the XML atom name
                pdb_index = atom['index']
                xml_atom_to_pdb_index[atom_name] = pdb_index
            
            # Generate CONECT records for bonds in this residue instance
            bonds_in_residue = 0
            for atom1_name, atom2_name in self.residue_bonds[residue_name]:
                if atom1_name in xml_atom_to_pdb_index and atom2_name in xml_atom_to_pdb_index:
                    idx1 = xml_atom_to_pdb_index[atom1_name]
                    idx2 = xml_atom_to_pdb_index[atom2_name]
                    
                    # Avoid duplicate bonds globally
                    bond_key = tuple(sorted([idx1, idx2]))
                    if bond_key not in bonds_added:
                        conect_records.append(f"CONECT{idx1:5}{idx2:5}\n")
                        bonds_added.add(bond_key)
                        bonds_in_residue += 1
                        total_bonds_processed += 1
            
            if bonds_in_residue > 0:
                print(f"Added {bonds_in_residue} bonds for {residue_name} (chain {chain}, seq {seq})")
        
        print(f"Total CONECT records generated: {total_bonds_processed}")
        return sorted(conect_records)


def main():
    """Main function to parse arguments and run the preprocessor."""
    parser = argparse.ArgumentParser(
        description="Preprocess PDB files for OpenMM compatibility"
    )
    parser.add_argument(
        "-pdb", 
        required=True, 
        help="Input PDB file path"
    )
    parser.add_argument(
        "-xml", 
        required=True, 
        help="Force field XML file path"
    )
    parser.add_argument(
        "-output", 
        required=True, 
        help="Output PDB file path"
    )
    
    args = parser.parse_args()
    
    # Validate input files
    if not Path(args.pdb).exists():
        print(f"Error: PDB file not found: {args.pdb}")
        sys.exit(1)
    
    if not Path(args.xml).exists():
        print(f"Error: XML file not found: {args.xml}")
        sys.exit(1)
    
    # Process the PDB file
    preprocessor = PDBPreprocessor(args.xml)
    preprocessor.process_pdb(args.pdb, args.output)


if __name__ == "__main__":
    main()