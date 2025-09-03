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

## What We Have Built

**✅ COMPLETED: Core Tool #1**
`off_to_openmm.py` - Complete .off file to OpenMM .xml converter
   ✅ Direct .off file parsing (handles MOLTYP/MOL, ATOMS/ATOM variations)
   ✅ Multi-molecule support with selective inclusion (-molnames flag)
   ✅ Special atom handling (4*, NETF, TORQ, M, MMM)
   ✅ Complete force field generation:
      - AtomTypes with element extraction
      - Residues with unique atom naming and connectivity
      - NonbondedForce with charges (sigma=0.0, epsilon=0.0, scales=1.0)
      - HarmonicBondForce with unit-converted parameters
      - HarmonicAngleForce with degree-to-radian conversion
      - PeriodicTorsionForce with phase conversion
      - CustomNonbondedForce for EXP interactions (exponential repulsion)
      - CustomNonbondedForce for SRD interactions (short-range dispersion, multiple powers)
   ✅ External charge file integration (-charges flag)
   ✅ Proper ForceField XML wrapper structure

**🔄 IN PROGRESS: Core Tool #2**
`pdb_preprocessor.py` - PDB file CONECT record generator
   - Adds bonding information to PDB files for OpenMM compatibility
   - Uses force field .xml file to determine connectivity

**✅ COMPLETED: Testing Infrastructure**
✅ Integration tests with real molecular systems
✅ Multiple test .off files: 1-butanol, hydrated_cyclohexene 
✅ Charge file testing with AtomType-Charge format
✅ pytest framework with unit and integration test structure

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
