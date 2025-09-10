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
        # Common residue name mappings
        self.residue_mapping = {
            'CYC': 'CYCQM',
            'SOL': 'H2OQM',
            'WAT': 'H2OQM'
        }
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
        
        # Process each line
        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                # Parse ATOM/HETATM record (PDB format specification)
                record_type = line[0:6].strip()
                atom_index = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()  # Standard PDB: columns 18-20 for residue name
                chain_id = line[21:22]
                residue_seq = line[22:26].strip()
                
                # Debug: print parsing info
                print(f"Debug: Parsed residue_name='{residue_name}' from line: {line.strip()}")
                
                # Check if atom name is a class and convert if needed
                original_atom_name = atom_name
                if self._is_atom_class(atom_name):
                    atom_name = self._convert_atom_name(atom_name)
                    # Update the line with new atom name
                    line = line[:12] + f"{atom_name:>4}" + line[16:]
                
                atom_records.append({
                    'line': line,
                    'index': atom_index,
                    'name': atom_name,
                    'original_name': original_atom_name,
                    'residue': residue_name,
                    'chain': chain_id,
                    'seq': residue_seq
                })
                
                # Map for CONECT generation
                key = (residue_name, atom_name)
                atom_index_map[key] = atom_index
                
            else:
                # Keep other records as-is
                other_records.append(line)
        
        # Generate CONECT records
        conect_records = self._generate_conect_records(atom_records, atom_index_map)
        
        # Write output PDB
        try:
            with open(output_pdb, 'w') as f:
                # Write ATOM/HETATM records
                for atom in atom_records:
                    f.write(atom['line'])
                
                # Write other records (excluding END if present)
                for record in other_records:
                    if not record.startswith('END'):
                        f.write(record)
                
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
        
        # Group atoms by residue
        residues = {}
        for atom in atom_records:
            res_key = (atom['residue'], atom['chain'], atom['seq'])
            if res_key not in residues:
                residues[res_key] = []
            residues[res_key].append(atom)
        
        # Generate CONECT records for each residue instance
        for (residue_name, chain, seq), atoms in residues.items():
            if residue_name not in self.residue_bonds:
                print(f"Warning: No bonds found for residue {residue_name}")
                continue
            
            # Create mapping from XML atom names to PDB indices
            # We need to map PDB atom types to XML atom names
            xml_atom_to_pdb_index = {}
            
            if residue_name in self.type_to_atoms:
                # Group PDB atoms by type
                pdb_atoms_by_type = {}
                for atom in atoms:
                    atom_type = atom['name']  # PDB atom name is actually the type
                    if atom_type not in pdb_atoms_by_type:
                        pdb_atoms_by_type[atom_type] = []
                    pdb_atoms_by_type[atom_type].append(atom)
                
                # Map XML atom names to PDB indices using type matching
                for atom_type, xml_atom_names in self.type_to_atoms[residue_name].items():
                    if atom_type in pdb_atoms_by_type:
                        pdb_atoms = pdb_atoms_by_type[atom_type]
                        # Map each XML atom name to a PDB atom of the same type
                        for i, xml_atom_name in enumerate(xml_atom_names):
                            if i < len(pdb_atoms):
                                xml_atom_to_pdb_index[xml_atom_name] = pdb_atoms[i]['index']
            
            # Generate CONECT records for bonds in this residue
            bonds_added = set()
            bonds_processed = 0
            for atom1_name, atom2_name in self.residue_bonds[residue_name]:
                if atom1_name in xml_atom_to_pdb_index and atom2_name in xml_atom_to_pdb_index:
                    idx1 = xml_atom_to_pdb_index[atom1_name]
                    idx2 = xml_atom_to_pdb_index[atom2_name]
                    
                    # Avoid duplicate bonds
                    bond_key = tuple(sorted([idx1, idx2]))
                    if bond_key not in bonds_added:
                        conect_records.append(f"CONECT{idx1:5}{idx2:5}\n")
                        bonds_added.add(bond_key)
                        bonds_processed += 1
            
            print(f"Processed {bonds_processed} bonds for residue {residue_name}")
        
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