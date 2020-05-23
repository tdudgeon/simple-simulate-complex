from simtk.openmm.app import *
from simtk.openmm import *
from simtk import unit
from sys import stdout
from pdbfixer import PDBFixer
from openmmforcefields.generators import SystemGenerator

fixer = PDBFixer(filename='protein_orig.pdb')
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.findNonstandardResidues()
print('Residues:', fixer.missingResidues)
print('Atoms:', fixer.missingAtoms)
print('Terminals:', fixer.missingTerminals)
print('Non-standard:', fixer.nonstandardResidues)

fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.4)
fixer.removeHeterogens(False)

with open('protein_prepared.pdb', 'w') as outfile:
    PDBFile.writeFile(fixer.topology, fixer.positions, file=outfile, keepIds=True)

system_generator = SystemGenerator(forcefields=['amber/ff14SB.xml'])
system = system_generator.create_system(fixer.topology)
integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
simulation = Simulation(fixer.topology, system, integrator)
simulation.context.setPositions(fixer.positions)
print('Minimising')
simulation.minimizeEnergy()

# write out the minimised PDB
with open('minimised1.pdb', 'w') as outfile:
    PDBFile.writeFile(fixer.topology, simulation.context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(), file=outfile, keepIds=True)

print('Done')