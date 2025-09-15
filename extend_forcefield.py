#!/usr/bin/env python3
"""
Force Field Extension Tool

Extends an existing OpenMM XML force field by adding new atoms, bonds, angles,
and dihedrals to create extended molecules (e.g., extending 1-butanol to n-octanol).

This tool preserves all atom types and nonbonded interactions from the base force field
while adding new bonded interactions for the extended molecule structure.
"""

import argparse
import xml.etree.ElementTree as ET
import sys
from collections import defaultdict


class ForceFieldExtender:
    """Class to handle extension of OpenMM XML force fields."""

    def __init__(self):
        self.base_ff = None
        self.new_atoms = []
        self.new_bonds = []
        self.extra_terms = {'angles': [], 'dihedrals': []}
        self.atom_type_mapping = {}

    def parse_atoms_file(self, atoms_file):
        """Parse atoms file containing atom indices and types."""
        print(f"Parsing atoms file: {atoms_file}")
        with open(atoms_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    parts = line.split()
                    if len(parts) >= 2:
                        atom_idx = int(parts[0])
                        atom_type = parts[1]
                        self.new_atoms.append((atom_idx, atom_type))
                        self.atom_type_mapping[atom_idx] = atom_type
        print(f"Loaded {len(self.new_atoms)} atoms")

    def parse_bonds_file(self, bonds_file):
        """Parse bonds file containing atom connectivity."""
        print(f"Parsing bonds file: {bonds_file}")
        with open(bonds_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    parts = line.split()
                    if len(parts) >= 2:
                        atom1 = int(parts[0])
                        atom2 = int(parts[1])
                        self.new_bonds.append((atom1, atom2))
        print(f"Loaded {len(self.new_bonds)} bonds")

    def parse_extra_terms_file(self, extra_terms_file):
        """Parse extra terms file containing angle and dihedral definitions."""
        print(f"Parsing extra terms file: {extra_terms_file}")
        with open(extra_terms_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#define'):
                    self._parse_define_line(line)
        print(f"Loaded {len(self.extra_terms['angles'])} angles and {len(self.extra_terms['dihedrals'])} dihedrals")

    def _parse_define_line(self, line):
        """Parse a #define line for bonded interactions."""
        parts = line.split()
        if len(parts) < 4:
            return

        atom_types = parts[1].split('_')
        interaction_type = parts[2]

        if interaction_type == 'HAR' and len(atom_types) == 3:
            # Angle: HAR angle_degrees force_constant
            angle_deg = float(parts[3])
            force_const = float(parts[4])
            angle_rad = angle_deg * 3.14159265359 / 180.0  # Convert to radians
            force_const_openmm = force_const * 4.184  # Convert kcal/mol to kJ/mol

            self.extra_terms['angles'].append({
                'types': atom_types,
                'angle': angle_rad,
                'k': force_const_openmm
            })

        elif interaction_type == 'NCO' and len(atom_types) == 4:
            # Dihedral: NCO phase_deg force_constant periodicity
            phase_deg = float(parts[3])
            force_const = float(parts[4])
            periodicity = int(parts[5])
            phase_rad = phase_deg * 3.14159265359 / 180.0  # Convert to radians
            # Use original force constant units from file

            self.extra_terms['dihedrals'].append({
                'types': atom_types,
                'phase': phase_rad,
                'k': force_const,
                'periodicity': periodicity
            })

    def load_base_forcefield(self, xml_file):
        """Load the base XML force field."""
        print(f"Loading base force field: {xml_file}")
        tree = ET.parse(xml_file)
        self.base_ff = tree.getroot()
        print("Base force field loaded successfully")

    def extend_residue(self):
        """Replace the residue definition with the new extended molecule structure."""
        residues = self.base_ff.find('Residues')
        if residues is None:
            raise ValueError("No Residues section found in base force field")

        # Find the first (and presumably only) residue
        residue = residues.find('Residue')
        if residue is None:
            raise ValueError("No Residue found in Residues section")

        # Clear existing atoms and bonds, keep residue name
        residue_name = residue.get('name', 'UNK')
        residues.remove(residue)

        # Create new residue with extended structure
        new_residue = ET.SubElement(residues, 'Residue')
        new_residue.set('name', residue_name)
        new_residue.text = '\n'

        # Create atom name mapping based on atom indices and types
        atom_name_map = {}
        carbon_count = 0
        hydrogen_count = 0
        oxygen_count = 0

        # First pass: Add all atoms (before bonds)
        for atom_idx, atom_type in self.new_atoms:
            if atom_type.startswith('C'):
                atom_name = f"C{carbon_count}"
                carbon_count += 1
            elif atom_type.startswith('H'):
                atom_name = f"H{hydrogen_count}"
                hydrogen_count += 1
            elif atom_type.startswith('O'):
                atom_name = f"O{oxygen_count}"
                oxygen_count += 1
            else:
                atom_name = f"X{atom_idx-1}"  # Fallback for unknown types

            atom_name_map[atom_idx] = atom_name

            # Add atom to residue
            atom_elem = ET.SubElement(new_residue, 'Atom')
            atom_elem.set('name', atom_name)
            atom_elem.set('type', atom_type)
            atom_elem.tail = '\n'

        # Second pass: Add all bonds (after atoms)
        for atom1_idx, atom2_idx in self.new_bonds:
            if atom1_idx in atom_name_map and atom2_idx in atom_name_map:
                bond_elem = ET.SubElement(new_residue, 'Bond')
                bond_elem.set('atomName1', atom_name_map[atom1_idx])
                bond_elem.set('atomName2', atom_name_map[atom2_idx])
                bond_elem.tail = '\n'

        print(f"Residue replaced with extended structure: {len(self.new_atoms)} atoms, {len(self.new_bonds)} bonds")

    def _get_atom_name_by_index(self, residue, atom_idx):
        """Get atom name by its index (1-based) in the residue."""
        atoms = residue.findall('Atom')
        if 1 <= atom_idx <= len(atoms):
            return atoms[atom_idx - 1].get('name')
        return None

    def add_bonded_interactions(self):
        """Add new angle and dihedral interactions to force field."""
        self._add_angles()
        self._add_dihedrals()

    def _add_angles(self):
        """Add new angle interactions."""
        angle_force = self.base_ff.find('HarmonicAngleForce')
        if angle_force is None:
            return

        for angle_term in self.extra_terms['angles']:
            # Check if this angle type already exists
            existing = False
            for existing_angle in angle_force.findall('Angle'):
                if (existing_angle.get('class1') == angle_term['types'][0] and
                    existing_angle.get('class2') == angle_term['types'][1] and
                    existing_angle.get('class3') == angle_term['types'][2]):
                    existing = True
                    break

            if not existing:
                angle_elem = ET.SubElement(angle_force, 'Angle')
                angle_elem.set('class1', angle_term['types'][0])
                angle_elem.set('class2', angle_term['types'][1])
                angle_elem.set('class3', angle_term['types'][2])
                angle_elem.set('angle', str(angle_term['angle']))
                angle_elem.set('k', str(angle_term['k']))

        print(f"Added {len(self.extra_terms['angles'])} new angle interactions")

    def _add_dihedrals(self):
        """Add new dihedral interactions."""
        dihedral_force = self.base_ff.find('PeriodicTorsionForce')
        if dihedral_force is None:
            return

        for dihedral_term in self.extra_terms['dihedrals']:
            # Check if this dihedral type already exists
            existing = False
            for existing_dihedral in dihedral_force.findall('Proper'):
                if (existing_dihedral.get('class1') == dihedral_term['types'][0] and
                    existing_dihedral.get('class2') == dihedral_term['types'][1] and
                    existing_dihedral.get('class3') == dihedral_term['types'][2] and
                    existing_dihedral.get('class4') == dihedral_term['types'][3]):
                    existing = True
                    break

            if not existing:
                dihedral_elem = ET.SubElement(dihedral_force, 'Proper')
                dihedral_elem.set('class1', dihedral_term['types'][0])
                dihedral_elem.set('class2', dihedral_term['types'][1])
                dihedral_elem.set('class3', dihedral_term['types'][2])
                dihedral_elem.set('class4', dihedral_term['types'][3])
                dihedral_elem.set('periodicity1', str(dihedral_term['periodicity']))
                dihedral_elem.set('phase1', str(dihedral_term['phase']))
                dihedral_elem.set('k1', str(dihedral_term['k']))
                dihedral_elem.tail = '\n'

        print(f"Added {len(self.extra_terms['dihedrals'])} new dihedral interactions")

    def write_extended_forcefield(self, output_file):
        """Write the extended force field to XML file."""
        # Fix formatting for all major sections
        self._format_xml_tree()

        # Create a clean XML tree with proper formatting
        tree = ET.ElementTree(self.base_ff)

        # Write to file with XML declaration
        with open(output_file, 'wb') as f:
            tree.write(f, encoding='utf-8', xml_declaration=True)

        print(f"Extended force field written to: {output_file}")

    def _format_xml_tree(self):
        """Add proper formatting to the XML tree."""
        # Add newlines after major closing tags
        for section in ['AtomTypes', 'Residues', 'NonbondedForce', 'HarmonicBondForce',
                       'HarmonicAngleForce', 'PeriodicTorsionForce', 'CustomNonbondedForce']:
            element = self.base_ff.find(section)
            if element is not None:
                element.tail = '\n\n'

        # Fix the final closing tag to have a newline
        last_element = list(self.base_ff)[-1]
        if last_element is not None:
            last_element.tail = '\n'


def main():
    parser = argparse.ArgumentParser(
        description="Extend an OpenMM XML force field with new atoms and bonded interactions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  ./extend_forcefield.py -atoms octanol_atoms.dat -bonds octanol_bonds.dat \\
                        -extend_ff 1-butanol.xml -extra_terms octanol_extra_terms.dat \\
                        -o octanol.xml
        """
    )

    parser.add_argument('-atoms', required=True,
                       help='File containing atom indices and types')
    parser.add_argument('-bonds', required=True,
                       help='File containing bond connectivity')
    parser.add_argument('-extend_ff', required=True,
                       help='Base XML force field to extend')
    parser.add_argument('-extra_terms', required=True,
                       help='File containing extra bonded terms (angles, dihedrals)')
    parser.add_argument('-o', '--output', default='extended_forcefield.xml',
                       help='Output XML file (default: extended_forcefield.xml)')

    args = parser.parse_args()

    try:
        # Initialize the force field extender
        extender = ForceFieldExtender()

        # Load all input data
        extender.parse_atoms_file(args.atoms)
        extender.parse_bonds_file(args.bonds)
        extender.parse_extra_terms_file(args.extra_terms)
        extender.load_base_forcefield(args.extend_ff)

        # Perform the extension
        extender.extend_residue()
        extender.add_bonded_interactions()

        # Write output
        extender.write_extended_forcefield(args.output)

        print("Force field extension completed successfully!")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()