import sys
from openforcefield.topology import Molecule
from openmmforcefields.generators import SystemGenerator
from simtk import unit
from simtk.openmm import app
from simtk.openmm.app import PDBFile, Simulation, Modeller, PDBReporter, StateDataReporter, DCDReporter
from simtk.openmm import *
import numpy as np

# check whether we have a GPU platform and if so set the precision to mixed
speed = 0
for i in range(Platform.getNumPlatforms()):
    p = Platform.getPlatform(i)
    # print(p.getName(), p.getSpeed())
    if p.getSpeed() > speed:
        platform = p
        speed = p.getSpeed()

if platform.getName() == 'CUDA' or platform.getName() == 'OpenCL':
    platform.setPropertyDefaultValue('Precision', 'mixed')
    print('Set precision for platform', platform.getName(), 'to mixed')


# Create an openforcefield Molecule object
ligand_mol = Molecule.from_file('ligand1.sdf', file_format='sdf')
print(ligand_mol)

print('Preparing system')
# Initialize a SystemGenerator using the GAFF for the ligand and tip3p for the water.
forcefield_kwargs = {'constraints': app.HBonds, 'rigidWater': True, 'removeCMMotion': False, 'hydrogenMass': 4*unit.amu }
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

# Solvate
print('Adding solvent...')
# we use the 'padding' option to define the periodic box. The PDB file does not contain any
# unit cell information so we just create a box that has a 10A padding around the complex.
modeller.addSolvent(system_generator.forcefield, model='tip3p', padding=10.0*unit.angstroms)
print('System has %d atoms' % modeller.topology.getNumAtoms())

with open('complex1.pdb', 'w') as outfile:
    PDBFile.writeFile(modeller.topology, modeller.positions, outfile)

# Create the system using the SystemGenerator
system = system_generator.create_system(modeller.topology, molecules=ligand_mol)
integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
print('Uses Periodic box:', system.usesPeriodicBoundaryConditions())
print('Default Periodic box:', system.getDefaultPeriodicBoxVectors())

simulation = Simulation(modeller.topology, system, integrator, platform=platform)
context = simulation.context
context.setPositions(modeller.positions)

print('Minimising')
simulation.minimizeEnergy()

# Write out the minimised PDB. The 'enforcePeriodicBox=False' bit is important otherwise the different
# components can end up in different periodic boxes resulting in really strange looking output.
with open('minimised1.pdb', 'w') as outfile:
    PDBFile.writeFile(modeller.topology, context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(), file=outfile, keepIds=True)

# run the simulation. The 3rd arg to PDBReporter is important. Again, this applies the
# 'enforcePeriodicBox=False' logic to ensure you get sensible output.
simulation.reporters.append(PDBReporter('output1.pdb', 1000, False))
simulation.reporters.append(StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, temperature=True))
print('Starting simulation')
simulation.step(500000)

print('Done')
