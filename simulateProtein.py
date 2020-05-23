from openmmforcefields.generators import SystemGenerator
from simtk import unit
from simtk.openmm.app import PDBFile, Simulation
from simtk.openmm import *

protein_pdb = PDBFile('protein.pdb')

print('Preparing system')
# Initialize a SystemGenerator
system_generator = SystemGenerator(forcefields=['amber/ff14SB.xml'])

system = system_generator.create_system(protein_pdb.topology)
integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
simulation = Simulation(protein_pdb.topology, system, integrator)
simulation.context.setPositions(protein_pdb.positions)
print('Minimising')
simulation.minimizeEnergy()

# write out the minimised PDB
with open('minimised1.pdb', 'w') as outfile:
    PDBFile.writeFile(protein_pdb.topology, simulation.context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(), file=outfile, keepIds=True)

print('Done')
