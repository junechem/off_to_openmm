# Development Tasks

## Completed 
- [x] Set up GitHub repository (private)
- [x] Create project documentation files
- [x] Set up testing infrastructure with pytest
- [x] Create test directory structure with sample .off files
- [x] Set up development tools (black, ruff)
- [x] Create empty script templates with command line interfaces
- [x] Analyze old script approaches and .off file structure
- [x] Update documentation with multi-molecule support requirements

## Current Priority =
- [ ] **Implement off_to_openmm.py with direct .off file parsing**
  - Parse .off files directly (not intermediate .dat files)
  - Handle multiple molecule types (UNK, CYCQM, etc.) 
  - Support format variations ([MOLTYP]/[MOL], [ATOMS]/[ATOM])
  - Extract bonded interactions (bonds, angles, dihedrals)
  - Extract nonbonded interactions (COU, EXP, SRD)
  - Parse #define parameter statements
  - Generate complete OpenMM XML with all force sections

## Pending Tasks =Ë
- [ ] Implement -molnames flag for selective molecule inclusion
- [ ] Implement -charges flag for external charge file integration  
- [ ] Add comprehensive error handling and validation
- [ ] Create integration tests with real molecular systems
- [ ] Implement pdb_preprocessor.py for CONECT record generation
- [ ] Add unit tests for individual parsing components
- [ ] Performance optimization for large .off files
- [ ] Documentation and usage examples

## Test Cases Available >ê
- 1-butanol.off + charges file
- hydrated_cyclohexene.off + charges file  
- Original intra.off (multi-molecule system)

## Technical Notes =Ý
- .off files contain complete force field definitions
- Multiple format variations need support
- Charges come from separate .dat files with "AtomType Charge" format
- Final parameters defined via #define statements at file end
- Nonbonded forces need CustomNonbondedForce XML sections