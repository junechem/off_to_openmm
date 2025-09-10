#!/usr/bin/env python3

import os
import sys
from collections import defaultdict

def check_required_files():
    """Check for required files and provide helpful error messages"""
    required_files = {
        'atoms.dat': 'Contains atom IDs and types (format: "1 O0")',
        'parameters.dat': 'Contains interaction parameters for SRD and EXP'
    }
    
    missing_files = []
    for filename, description in required_files.items():
        if not os.path.exists(filename):
            missing_files.append((filename, description))
    
    if missing_files:
        print("ERROR: Required files are missing!")
        print("\nThe following files are required in the current directory:")
        for filename, description in missing_files:
            print(f"  - {filename}: {description}")
        
        print("\nOptional files:")
        print("  - charges.dat: Contains atomic charges (format: \"C1 0.24499\")")
        print("    Note: If missing or all charges are zero, charge section will be skipped")
        
        return False
    
    return True

# ============================================================================
# SHARED FUNCTION - read_atoms
# ============================================================================

def read_atoms(filename):
    """Read atoms.dat file and return unique atom types with indices"""
    unique_types = []
    
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                atom_type = parts[1]
                if atom_type not in unique_types:
                    unique_types.append(atom_type)
    
    return unique_types

# ============================================================================
# CHARGES SECTION - From generate_charges_xml.py
# ============================================================================

def read_charges(filename):
    """Read charges.dat file and return charge mapping"""
    charges = {}
    
    if not os.path.exists(filename):
        print(f"Warning: {filename} not found. Skipping charge generation.")
        return None
    
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                atom_type = parts[0]
                charge = float(parts[1])
                charges[atom_type] = charge
    
    return charges

def check_all_charges_zero(charges):
    """Check if all charges are zero (within floating point tolerance)"""
    tolerance = 1e-10
    return all(abs(charge) < tolerance for charge in charges.values())

def generate_charges_xml(charges, atom_types):
    """Generate charge information (NonbondedForce section removed)"""
    
    # This function now returns empty string as NonbondedForce section is not needed
    return ""

# ============================================================================
# SRD SECTION - From generate_srd_xml.py
# ============================================================================

def parse_srd_parameters(filename):
    """Parse parameters.dat file and extract SRD parameters grouped by power"""
    srd_parameters = defaultdict(list)
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            parts = line.split()
            if len(parts) >= 6 and parts[2] == 'SRD':
                atom1 = parts[0]
                atom2 = parts[1]
                p1_value = float(parts[3])  # P1: kcal*angstrom^(-P2)
                power = float(parts[4])     # P2: power
                p3_value = float(parts[5])  # P3: Angstrom^(-1)
                
                # Store parameters grouped by power
                srd_parameters[power].append({
                    'atom1': atom1,
                    'atom2': atom2,
                    'p1': p1_value,
                    'power': power,
                    'p3': p3_value
                })
    
    return srd_parameters

def convert_units(p1_value, power, p3_value):
    """Convert units from kcal/Angstrom to kJ/nm
    
    Based on the example values, there appears to be an additional factor needed.
    From comparing the example -0.004303143 vs calculated -0.002165441752,
    there's a factor of approximately 1.987 ≈ 2 needed.
    """
    # Unit conversions
    kcal_to_kj = 4.184  # kcal to kJ
    angstrom_to_nm = 0.1  # Angstrom to nm
    
    # P1: kcal*angstrom^(-P2) to kJ*nm^(-P2)
    
    p1_converted = p1_value * kcal_to_kj * (angstrom_to_nm**(-power))
    
    # P3: Angstrom to nm  
    p3_converted = p3_value * angstrom_to_nm
    
    return p1_converted, p3_converted

def create_interaction_matrices(srd_params_by_power, atom_types):
    """Create symmetric matrices for dispersion and r0 values for each power"""
    matrices_by_power = {}
    n_types = len(atom_types)
    
    for power, params in srd_params_by_power.items():
        # Initialize matrices with zeros
        disp_matrix = [[0.0 for _ in range(n_types)] for _ in range(n_types)]
        r0_matrix = [[0.0 for _ in range(n_types)] for _ in range(n_types)]
        
        # Fill matrices with parameters
        for param in params:
            atom1 = param['atom1']
            atom2 = param['atom2']
            
            if atom1 in atom_types and atom2 in atom_types:
                i = atom_types.index(atom1)
                j = atom_types.index(atom2)
                
                # Convert units
                p1_converted, p3_converted = convert_units(param['p1'], param['power'], param['p3'])
                
                # Fill symmetric matrices
                disp_matrix[i][j] = p1_converted
                disp_matrix[j][i] = p1_converted  # Symmetric
                
                r0_matrix[i][j] = p3_converted
                r0_matrix[j][i] = p3_converted  # Symmetric
        
        matrices_by_power[power] = {
            'disp_matrix': disp_matrix,
            'r0_matrix': r0_matrix
        }
    
    return matrices_by_power

def matrix_to_xml_string(matrix):
    """Convert matrix to XML string format (flattened)"""
    result = []
    for row in matrix:
        result.extend([str(val) for val in row])
    return '\t'.join(result)

def generate_custom_nonbonded_force_xml(power, matrices, atom_types):
    """Generate CustomNonbondedForce XML section for a specific power"""
    n_types = len(atom_types)
    disp_matrix = matrices['disp_matrix']
    r0_matrix = matrices['r0_matrix']
    
    # Convert power to integer for display if it's a whole number
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

# ============================================================================
# EXP SECTION - From generate_exp_xml.py
# ============================================================================

def parse_exp_parameters(filename):
    """Parse parameters.dat file and extract EXP parameters"""
    exp_parameters = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            parts = line.split()
            if len(parts) >= 5 and parts[2] == 'EXP':
                atom1 = parts[0]
                atom2 = parts[1]
                a_value = float(parts[3])      # P1: A (kcal/mol)
                alpha_value = float(parts[4])  # P2: alpha (Angstrom^-1)
                
                exp_parameters.append({
                    'atom1': atom1,
                    'atom2': atom2,
                    'a': a_value,
                    'alpha': alpha_value
                })
    
    return exp_parameters

def convert_exp_units(a_value, alpha_value):
    """Convert EXP units from kcal/mol and Angstrom^-1 to kJ/mol and nm^-1"""
    # Unit conversions
    kcal_to_kj = 4.184  # kcal to kJ
    angstrom_to_nm = 0.1  # Angstrom to nm
    
    # A: kcal/mol to kJ/mol
    a_converted = a_value * kcal_to_kj
    
    # alpha: Angstrom^-1 to nm^-1
    # Angstrom^-1 * (angstrom/0.1 nm)^-1 = Angstrom^-1 * 10
    alpha_converted = alpha_value / angstrom_to_nm
    
    return a_converted, alpha_converted

def create_exp_matrices(exp_parameters, atom_types):
    """Create symmetric matrices for A and alpha values"""
    n_types = len(atom_types)
    
    # Initialize matrices with zeros
    a_matrix = [[0.0 for _ in range(n_types)] for _ in range(n_types)]
    alpha_matrix = [[0.0 for _ in range(n_types)] for _ in range(n_types)]
    
    # Fill matrices with parameters
    for param in exp_parameters:
        atom1 = param['atom1']
        atom2 = param['atom2']
        
        if atom1 in atom_types and atom2 in atom_types:
            i = atom_types.index(atom1)
            j = atom_types.index(atom2)
            
            # Convert units
            a_converted, alpha_converted = convert_exp_units(param['a'], param['alpha'])
            
            # Fill symmetric matrices
            a_matrix[i][j] = a_converted
            a_matrix[j][i] = a_converted  # Symmetric
            
            alpha_matrix[i][j] = alpha_converted
            alpha_matrix[j][i] = alpha_converted  # Symmetric
    
    return a_matrix, alpha_matrix

def generate_exp_xml(a_matrix, alpha_matrix, atom_types):
    """Generate CustomNonbondedForce XML section for EXP interactions"""
    n_types = len(atom_types)
    
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

# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main():
    """Main function to generate unified nonbonded XML"""
    
    print("OpenMM Nonbonded Force Field Generator")
    print("=" * 40)
    
    # Check required files
    if not check_required_files():
        sys.exit(1)
    
    # Read atom types
    print("\n1. Reading atom types...")
    atom_types = read_atoms('atoms.dat')
    print(f"   Found {len(atom_types)} unique atom types: {atom_types}")
    
    xml_sections = []
    
    # ========================================================================
    # CHARGES SECTION (removed - charges now handled in residue section)
    # ========================================================================
    print("\n2. Processing charges...")
    charges = read_charges('charges.dat')
    
    if charges is None:
        print("   charges.dat not found - charges will be handled in residue section")
    elif check_all_charges_zero(charges):
        print("   All charges are zero - charges will be handled in residue section")
    else:
        print("   NonbondedForce section removed - charges handled in residue section")
        
        # NonbondedForce section is no longer generated
        # charges_xml = generate_charges_xml(charges, atom_types)
        # xml_sections.append(charges_xml)
    
    # ========================================================================
    # SRD SECTION
    # ========================================================================
    print("\n3. Processing SRD interactions...")
    srd_params_by_power = parse_srd_parameters('parameters.dat')
    
    if not srd_params_by_power:
        print("   No SRD parameters found - skipping SRD section")
    else:
        powers = sorted(srd_params_by_power.keys())
        print(f"   Found SRD parameters for powers: {powers}")
        
        # Create interaction matrices for each power
        matrices_by_power = create_interaction_matrices(srd_params_by_power, atom_types)
        
        # Generate XML for each power
        for power in powers:
            print(f"   Generating CustomNonbondedForce for SRD power {power}")
            srd_xml = generate_custom_nonbonded_force_xml(power, matrices_by_power[power], atom_types)
            xml_sections.append(srd_xml)
    
    # ========================================================================
    # EXP SECTION
    # ========================================================================
    print("\n4. Processing EXP interactions...")
    exp_parameters = parse_exp_parameters('parameters.dat')
    
    if not exp_parameters:
        print("   No EXP parameters found - skipping EXP section")
    else:
        print(f"   Found {len(exp_parameters)} EXP parameter sets")
        
        # Create interaction matrices
        a_matrix, alpha_matrix = create_exp_matrices(exp_parameters, atom_types)
        
        # Generate XML
        print("   Generating CustomNonbondedForce for EXP interactions")
        exp_xml = generate_exp_xml(a_matrix, alpha_matrix, atom_types)
        xml_sections.append(exp_xml)
    
    # ========================================================================
    # WRITE OUTPUT
    # ========================================================================
    print("\n5. Writing output...")
    
    if not xml_sections:
        print("   No force sections generated - no output file created")
        print("   This could happen if:")
        print("   - charges.dat is missing or all charges are zero")
        print("   - parameters.dat contains no SRD or EXP entries")
    else:
        output_filename = 'nonbonded.xml'
        with open(output_filename, 'w') as f:
            f.write('\n\n'.join(xml_sections) + '\n')
        
        print(f"   Generated {output_filename} with {len(xml_sections)} force section(s)")
        print(f"   File size: {os.path.getsize(output_filename)} bytes")
    
    print("\nDone!")

if __name__ == "__main__":
    main()