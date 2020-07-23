import sys
from simtk.openmm.app import *
from simtk.openmm import *
from simtk import unit
from pdbfixer import PDBFixer
from openmmforcefields.generators import SystemGenerator

if len(sys.argv) != 3:
    print('Usage: python prepareProtein.py input.pdb output')
    print('Creates outputs named output_fixed.pdb and output_minimised.pdb')
    exit(1)

pdb_in = sys.argv[1]
pdb_out = sys.argv[2]
print('Processing', pdb_in, 'to', pdb_out)

fixer = PDBFixer(filename=pdb_in)
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

with open(pdb_out + '_fixed.pdb', 'w') as outfile:
    PDBFile.writeFile(fixer.topology, fixer.positions, file=outfile, keepIds=True)

system_generator = SystemGenerator(forcefields=['amber/ff14SB.xml'])
system = system_generator.create_system(fixer.topology)
integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
simulation = Simulation(fixer.topology, system, integrator)
simulation.context.setPositions(fixer.positions)
print('Minimising')
simulation.minimizeEnergy()

# write out the minimised PDB
with open(pdb_out + '_minimised.pdb', 'w') as outfile:
    PDBFile.writeFile(fixer.topology, simulation.context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(), file=outfile, keepIds=True)

print('Done')