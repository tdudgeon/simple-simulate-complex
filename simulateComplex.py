import sys
from openforcefield.topology import Molecule
from openmmforcefields.generators import SystemGenerator
from simtk import unit
from simtk.openmm import app
from simtk.openmm.app import PDBFile, Simulation, Modeller
from simtk.openmm import *
import parmed

# Create an openforcefield Molecule object
ligand_mol = Molecule(open('ligand1.mol', 'rb'), file_format='mol')
# can't read as PDB as "No toolkits in registry can read file"
#complex_pdb = Molecule(open('complex1.pdb', 'rb'), file_format='pdb')

# Use Modeller to combine the protein and ligand into a complex
protein_pdb = PDBFile('protein.pdb')
ligand_pdb = PDBFile('ligand1.pdb')

# Approach 1
modeller = Modeller(protein_pdb.topology, protein_pdb.positions)
modeller.add(ligand_pdb.topology, ligand_pdb.positions)
PDBFile.writeFile(modeller.topology, modeller.positions, open('complex1.pdb', 'w'))

# Approach 2
# protein_pdb = parmed.load_file('protein.pdb')
# ligand_pdb = parmed.load_file('ligand1.pdb')
# complex = protein_pdb + ligand_pdb
# complex.save('complex1.pdb')

complex_pdb = PDBFile('complex1.pdb')

print('Preparing system')
forcefield_kwargs = { 'constraints' : app.HBonds, 'rigidWater' : True, 'removeCMMotion' : False, 'hydrogenMass' : 4*unit.amu }
# Initialize a SystemGenerator using GAFF
system_generator = SystemGenerator(
    forcefields=['amber/ff14SB.xml', 'amber/tip3p_standard.xml'],
    small_molecule_forcefield='gaff-2.11',
    forcefield_kwargs=forcefield_kwargs)

system = system_generator.create_system(complex_pdb.topology, molecules=[ligand_mol])
integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
simulation = Simulation(complex_pdb.topology, system, integrator)
simulation.context.setPositions(complex.positions)
print('Minimising')
simulation.minimizeEnergy()

print('Done')

