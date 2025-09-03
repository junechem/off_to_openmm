# Project Requirements Document (PRD)

## Problem Statement

**Primary Challenge:**
Currently, there exists no framework for converting force field parameters from CRYOFF (Create Your Own Force Field) output (.off files) to OpenMM .xml force field files. Researchers using AFM (Adaptive Force Matching) must manually copy all parameters to .xml files by hand, which is extremely tedious, time-consuming, and error-prone.

**Secondary Challenge:**
Standard .pdb files lack the bonding information (CONECT records) required by OpenMM for simulations with custom forces, necessitating manual editing of PDB files.

**Scope Expansion:**
The solution must handle multiple diverse molecular systems, including:
- Different molecule types (UNK, CYCQM, etc.)
- Various .off file format variations
- Different molecular structures (1-butanol, hydrated_cyclohexene, etc.)
- Separate charge assignment files for each system

## What We Are Building

**Core Tools:**
1. `off_to_openmm.py` - Direct .off file to OpenMM .xml converter
   - Parses .off files directly (not intermediate .dat files)
   - Supports multiple molecule types in single .off file
   - Handles bonded interactions (bonds, angles, dihedrals)
   - Handles nonbonded interactions (Coulombic, exponential, dispersion)
   - Allows selective molecule inclusion via command line flags

2. `pdb_preprocessor.py` - PDB file CONECT record generator
   - Adds bonding information to PDB files for OpenMM compatibility
   - Uses force field .xml file to determine connectivity

**Testing Infrastructure:**
- Integration tests with real molecular systems
- Multiple test .off files representing different chemical systems
- Charge file testing with various molecular configurations

## Target Audience

**Primary Users:**
- Computational chemists and physicists using AFM-generated force fields
- Researchers requiring GPU-accelerated OpenMM simulations
- Scientists performing Path Integral Molecular Dynamics (PIMD) with nuclear quantum effects

**Molecular Systems:**
- Organic molecules (alcohols, hydrocarbons)
- Solvated systems (molecules with explicit water)
- Custom force field applications requiring high precision

## Strategic Importance

This tool suite eliminates the current bottleneck preventing AFM force field adoption in OpenMM, enabling:
- Faster research workflows for computational chemistry
- Broader adoption of AFM methodology
- GPU-accelerated simulations with custom force fields
- Advanced simulation techniques (PIMD) with AFM models
