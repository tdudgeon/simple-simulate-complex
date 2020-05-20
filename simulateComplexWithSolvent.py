# This simulation is not working correctly. The minimisation step blows up.
# Probably something related to the periodic box boundaries.

import sys
from openforcefield.topology import Molecule
from openmmforcefields.generators import SystemGenerator
from simtk import unit
from simtk.openmm import app
from simtk.openmm.app import PDBFile, Simulation, Modeller, PDBReporter, StateDataReporter, DCDReporter
from simtk.openmm import *
import numpy as np

# Create an openforcefield Molecule object
ligand_mol = Molecule.from_file('ligand1.sdf', file_format='sdf')
print(ligand_mol)

print('Preparing system')
forcefield_kwargs = { 'constraints' : app.HBonds, 'rigidWater' : True, 'removeCMMotion' : False, 'hydrogenMass' : 4*unit.amu }
# Initialize a SystemGenerator using GAFF
system_generator = SystemGenerator(
    forcefields=['amber/ff14SB.xml', 'amber/tip3p_standard.xml'],
    small_molecule_forcefield='gaff-2.11',
    molecules=[ligand_mol],
    forcefield_kwargs=forcefield_kwargs)

# Use Modeller to combine the protein and ligand into a complex
print('Reading protein')
protein_pdb = PDBFile('protein.pdb')

# reading the ligand gives lots of warnings about "duplicate atom" but this seems incorrect
print('Reading ligand')
ligand_pdb = PDBFile('ligand1.pdb')

print('Preparing complex')
modeller = Modeller(protein_pdb.topology, protein_pdb.positions)
print('System has %d atoms' % modeller.topology.getNumAtoms())
modeller.add(ligand_pdb.topology, ligand_pdb.positions)
print('System has %d atoms' % modeller.topology.getNumAtoms())

box_vectors = unit.Quantity(np.diag([70, 70, 70]), unit.angstrom)
modeller.topology.setPeriodicBoxVectors(box_vectors)

# Solvate
print('Adding solvent...')
modeller.addSolvent(system_generator.forcefield, model='tip3p')#, padding=5.0*unit.angstroms)
print('System has %d atoms' % modeller.topology.getNumAtoms())

PDBFile.writeFile(modeller.topology, modeller.positions, open('complex1.pdb', 'w'))

complex_pdb = PDBFile('complex1.pdb')

system = system_generator.create_system(complex_pdb.topology, molecules=ligand_mol)
integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
print('Uses Periodic box:', system.usesPeriodicBoundaryConditions())
print('Default Periodic box:', system.getDefaultPeriodicBoxVectors())

simulation = Simulation(complex_pdb.topology, system, integrator)
context = simulation.context
context.setPositions(complex_pdb.positions)

print('Minimising')
simulation.minimizeEnergy()

with open('minimised1.pdb', 'w') as outfile:
    PDBFile.writeFile(modeller.topology, context.getState(getPositions=True,enforcePeriodicBox=True).getPositions(), file=outfile, keepIds=True)


simulation.reporters.append(PDBReporter('output1.pdb', 1000))
simulation.reporters.append(StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, temperature=True))
print('Starting simulation')
simulation.step(10000)

print('Done')
