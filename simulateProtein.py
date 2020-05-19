from openmmforcefields.generators import SystemGenerator
from simtk import unit
from simtk.openmm import app
from simtk.openmm.app import PDBFile, Simulation
from simtk.openmm import *

protein_pdb = PDBFile('protein.pdb')

print('Preparing system')
forcefield_kwargs = { 'constraints' : app.HBonds, 'rigidWater' : True, 'removeCMMotion' : False, 'hydrogenMass' : 4*unit.amu }
# Initialize a SystemGenerator using GAFF
system_generator = SystemGenerator(
    forcefields=['amber/ff14SB.xml', 'amber/tip3p_standard.xml'],
    small_molecule_forcefield='gaff-2.11',
    forcefield_kwargs=forcefield_kwargs)

system = system_generator.create_system(protein_pdb.topology)
integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
simulation = Simulation(protein_pdb.topology, system, integrator)
simulation.context.setPositions(protein_pdb.positions)
print('Minimising')
simulation.minimizeEnergy()

print('Done')

