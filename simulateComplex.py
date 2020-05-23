import sys
from openforcefield.topology import Molecule
from openmmforcefields.generators import SystemGenerator
from simtk import unit
from simtk.openmm import app
from simtk.openmm.app import PDBFile, Simulation, Modeller, PDBReporter, StateDataReporter
from simtk.openmm import *
import parmed

# Create an openforcefield Molecule object
ligand_mol = Molecule.from_file('ligand1.sdf', file_format='sdf')
print(ligand_mol)
# can't read as PDB as "No toolkits in registry can read file"
#complex_pdb = Molecule(open('complex1.pdb', 'rb'), file_format='pdb')

print('Reading protein')
protein_pdb = PDBFile('protein.pdb')

# reading the ligand gives lots of warnings about "duplicate atom" but this seems incorrect
print('Reading ligand')
ligand_pdb = PDBFile('ligand1.pdb')

# Use Modeller to combine the protein and ligand into a complex
print('Preparing complex')
modeller = Modeller(protein_pdb.topology, protein_pdb.positions)
modeller.add(ligand_pdb.topology, ligand_pdb.positions)
print('System has %d atoms' % modeller.topology.getNumAtoms())
PDBFile.writeFile(modeller.topology, modeller.positions, open('complex1.pdb', 'w'))

complex_pdb = PDBFile('complex1.pdb')

# Initialize a SystemGenerator using the GAFF for the ligand
print('Preparing system')
forcefield_kwargs = { 'constraints': app.HBonds, 'rigidWater': True, 'removeCMMotion': False, 'hydrogenMass': 4*unit.amu }
system_generator = SystemGenerator(
    forcefields=['amber/ff14SB.xml'],
    small_molecule_forcefield='gaff-2.11',
    forcefield_kwargs=forcefield_kwargs)

system = system_generator.create_system(complex_pdb.topology, molecules=ligand_mol)
integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
simulation = Simulation(complex_pdb.topology, system, integrator)
simulation.context.setPositions(complex_pdb.positions)
print('Minimising')
simulation.minimizeEnergy()

# write out the minimised PDB
with open('minimised1.pdb', 'w') as outfile:
    PDBFile.writeFile(complex_pdb.topology, simulation.context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(), file=outfile, keepIds=True)

# Run the sumulation
simulation.reporters.append(PDBReporter('output1.pdb', 1000))
simulation.reporters.append(StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, temperature=True))
print('Starting simulation')
simulation.step(10000)

print('Done')
