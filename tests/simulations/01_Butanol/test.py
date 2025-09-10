#!/usr/bin/env python3
from math import pi
from openmm.app import *
from openmm import *
from openmm.unit import *
import openmm as mm
import numpy as np
import xml.etree.ElementTree as ET
import re


### Simulation Control:
temperature     = 298*kelvin
pressure        = 1*bar
length          = 10  # picoseconds
timestep_size   = 0.0005 # picoseconds
friction_coeff  = 2/picosecond
dispcorr_freq   = 25  # steps
nmol            = 125  # molecules
rcut            = 1E-9  # meters

### r^-6 Dispersion Correction Function
def get_r6_correction(forcefield_xml, topology, cutoff_nm):
    """Extract r^-6 dispersion parameters and return correction function"""
    tree = ET.parse(forcefield_xml)
    root = tree.getroot()
    
    # Find r^-6 force and extract dispTable
    for force in root.findall('.//CustomNonbondedForce'):
        if 'r^6' in force.get('energy', ''):
            disp_func = force.find('.//Function[@name="dispTable"]')
            values = [float(x) for x in disp_func.text.split()]
            n_types = int(disp_func.get('xsize'))
            disp_table = np.array(values).reshape(n_types, n_types)
            
            # Get atom type mapping
            type_to_index = {}
            for atom in force.findall('.//Atom'):
                atom_type = atom.get('type')
                index = int(atom.get('t'))
                type_to_index[atom_type] = index
            break
    
    # Get atom name to type mapping from residue definition
    name_to_type = {}
    for residue in root.findall('.//Residue'):
        for atom_def in residue.findall('.//Atom'):
            name = atom_def.get('name')
            atom_type = atom_def.get('type')
            name_to_type[name] = atom_type
    
    # Count atoms of each type in the system
    atom_counts = np.zeros(n_types)
    for atom in topology.atoms():
        # Map atom name to forcefield type
        atom_name = atom.name
        if atom_name in name_to_type:
            atom_type = name_to_type[atom_name]
            if atom_type in type_to_index:
                atom_counts[type_to_index[atom_type]] += 1
    
    # Try a simpler approach: calculate total C6 per molecule-molecule interaction
    # Each molecule has these atom counts per molecule
    mol_atom_counts = atom_counts / nmol  # atoms of each type per molecule
    
    print(f"Atoms per molecule: {dict(zip(['O0','C1','C2','H0','H1','H2'], mol_atom_counts))}")
    
    # Calculate total C6 interaction between two molecules
    total_c6_per_mol_pair = 0
    pair_count = 0
    
    for i in range(n_types):
        for j in range(n_types):
            if mol_atom_counts[i] > 0 and mol_atom_counts[j] > 0:
                c6_ij = abs(disp_table[i][j])
                if c6_ij > 1e-10:  # Only non-zero interactions
                    # Number of i-j pairs between two molecules
                    n_pairs = mol_atom_counts[i] * mol_atom_counts[j]
                    total_c6_per_mol_pair += c6_ij * n_pairs
                    pair_count += n_pairs
    
    # This gives the total C6 interaction strength between two molecules
    # For the dispersion correction, we might need to scale this appropriately
    avg_c6 = total_c6_per_mol_pair
    
    print(f"Total C6 per molecule pair: {total_c6_per_mol_pair:.6e}")
    print(f"Number of interacting atom pairs per mol-mol: {pair_count}")
    
    print(f"Debug: calculated avg_c6 = {avg_c6:.6e} kJ*nm^6/mol")
    
    # Use the same formula as the working water example
    # disp_pam = raw_c6 * 1000 * (1E-9)^6 to convert kJ*nm^6/mol to J*m^6/mol
    disp_pam = avg_c6 * 1000 * ((1e-9)**6)
    avogadro = 6.02214129e23
    cutoff_m = cutoff_nm * 1e-9
    
    print(f"Effective disp_pam = {disp_pam:.6e} J*m^6/mol")
    
    def pcorr(volume_m3):
        """Return pressure correction using same formula as water example"""
        dens = nmol / volume_m3  # molecules/m^3
        sc = (-4 * pi / 3)
        correction = sc * disp_pam * dens**2 * (cutoff_m**(-3)) / avogadro / 1e5
        #print(f"Debug: dens={dens:.6e}, correction={correction:.6f} bar")
        return correction
    
    
    return pcorr



pdb = PDBFile('conf.pdb')
forcefield = ForceField('forcefield.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, nonbondedCutoff=rcut*meter)


print("---\nConfiguring Custom Nonbonded Forces---")
for force in system.getForces():
    if isinstance(force, mm.CustomNonbondedForce):
        original_energy = force.getEnergyFunction()
        print(f"\nFound CustomNonbondedForce with energy: {original_energy[:40]}...")

        # Set the cutoff method and distance for all custom forces
        force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        force.setCutoffDistance(rcut*meter) # Use rcut*meter to match createSystem call
        print(f"  - Set cutoff distance.")

        # Enable long-range correction for dispersion forces, disable for others.
        if "disp/" in original_energy:
            print("  - This is a dispersion force. Enabling long-range correction.")
            force.setUseLongRangeCorrection(True)
        else:
            print("  - This is not a dispersion force. Disabling long-range correction.")
            force.setUseLongRangeCorrection(False)
print("\n--- Custom Force Configuration Complete ---")






# Add NonbondedForce for PME electrostatics (matching CustomNonbondedForce exclusions)
nonbonded = NonbondedForce()
nonbonded.setNonbondedMethod(NonbondedForce.PME)
nonbonded.setCutoffDistance(rcut*meter)

# Add exclusions matching CustomNonbondedForce bondCutoff="2" (up to 2 bonds away)
# Build bond graph efficiently
atom_list = list(pdb.topology.atoms())
atom_to_index = {atom: i for i, atom in enumerate(atom_list)}
bond_graph = {i: set() for i in range(len(atom_list))}

# Build adjacency list
for bond in pdb.topology.bonds():
    i, j = atom_to_index[bond.atom1], atom_to_index[bond.atom2]
    bond_graph[i].add(j)
    bond_graph[j].add(i)

# Find all pairs within 2 bonds
exclusions = set()
for i in range(len(atom_list)):
    # 1 bond away (directly bonded)
    for j in bond_graph[i]:
        if i < j:
            exclusions.add((i, j))
    
    # 2 bonds away (1-3 interactions)
    for j in bond_graph[i]:
        for k in bond_graph[j]:
            if i < k and k not in bond_graph[i]:
                exclusions.add((i, k))

# Add exceptions
for i, j in exclusions:
    nonbonded.addException(i, j, 0.0*elementary_charge**2, 0.0*nanometer, 0.0*kilojoule_per_mole)

system.addForce(nonbonded)
n_beads = 32
system.addForce(RPMDMonteCarloBarostat(pressure))

# Setup dispersion correction
pcorr = get_r6_correction('forcefield.xml', pdb.topology, rcut*1e9)

for f in system.getForces():
    print(type(f).__name__)



for i, f in enumerate(system.getForces()):
    f.setForceGroup(i)

integrator = RPMDIntegrator(n_beads, temperature, friction_coeff, timestep_size)
platform = Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

print("=== INITIAL ENERGIES ===")
for i, f in enumerate(system.getForces()):
    state = simulation.context.getState(getEnergy=True, groups={i})
    energy = state.getPotentialEnergy()
    if isinstance(f, NonbondedForce):
        print(f"COULOMBIC (PME): {energy}")
    else:
        print(f"{f.getName()}: {energy}")

simulation.context.setVelocitiesToTemperature(temperature)

print("\n=== INITIAL ENERGIES (after setting velocities) ===")
for i, f in enumerate(system.getForces()):
    state = simulation.context.getState(getEnergy=True, groups={i})
    energy = state.getPotentialEnergy()
    volume = state.getPeriodicBoxVolume()
    if isinstance(f, NonbondedForce):
        print(f"COULOMBIC (PME): {energy}")
    else:
        print(f"{f.getName()}: {energy}")
print(f"Initial Volume: {volume}")

simulation.reporters.append(StateDataReporter("md.log", 2000, time = True, potentialEnergy = True, kineticEnergy = True, totalEnergy = True, temperature = True, density=True, separator='\t', volume=True, speed=True))
simulation.reporters.append(DCDReporter('traj.dcd', 1000))
print('Using Platform:', simulation.context.getPlatform().getName())

for i in np.arange(0, length/timestep_size, dispcorr_freq):
    volume_m3 = simulation.context.getState().getPeriodicBoxVolume().value_in_unit(meter**3)
    pressure_correction = pcorr(volume_m3)
    new_p = 1 - pressure_correction  # Same pattern as working Water_Example
    #print(f"Volume: {volume_m3:.6e} m³, Pcorr: {pressure_correction:.6f} bar, New P: {new_p:.6f} bar")
    simulation.context.setParameter(RPMDMonteCarloBarostat.Pressure(), new_p*bar)
    simulation.step(dispcorr_freq)

print("\n=== FINAL ENERGIES (after simulation) ===")
for i, f in enumerate(system.getForces()):
    state = simulation.context.getState(getEnergy=True, groups={i})
    energy = state.getPotentialEnergy()
    if isinstance(f, NonbondedForce):
        print(f"COULOMBIC (PME): {energy}")
    else:
        print(f"{f.getName()}: {energy}")

# After simulation.step(...) is complete
positions = simulation.context.getState(getPositions=True).getPositions()

with open('confout.pdb', 'w') as f:
    PDBFile.writeFile(simulation.topology, positions, f)

