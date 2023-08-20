import sys, time
from openforcefield.topology import Molecule
from openmmforcefields.generators import SystemGenerator
from simtk import unit, openmm
from simtk.openmm import app, LangevinIntegrator
from simtk.openmm.app import PDBFile, Simulation, Modeller, PDBReporter, DCDReporter, StateDataReporter

from rdkit import Chem

import utils

t0 = time.time()

temperature = 330 * unit.kelvin
equilibration_steps = 200
reporting_interval = 1000

if len(sys.argv) != 5:
    print('Usage: python simulateComplexWithSolvent2.py input.pdb input.mol output num_steps')
    print('Prepares complex of input.pdb and input.mol and generates complex named output_complex.pdb,')
    print(' minimised complex named output_minimised.pdb and MD trajectory named output_traj.pdb and/or output_traj.dcd')
    exit(1)


pdb_in = sys.argv[1]
mol_in = sys.argv[2]
output_complex = sys.argv[3] + '_complex.pdb'
output_traj_pdb = sys.argv[3] + '_traj.pdb'
output_traj_dcd = sys.argv[3] + '_traj.dcd'
output_min = sys.argv[3] + '_minimised.pdb'
num_steps = int(sys.argv[4])
print('Processing', pdb_in, 'and', mol_in, 'with', num_steps, 'steps generating outputs',
      output_complex, output_min, output_traj_pdb, output_traj_dcd)

# get the chosen or fastest platform
platform = utils.get_platform()

# Read the molfile into RDKit, add Hs and create an openforcefield Molecule object
print('Reading ligand')
rdkitmol = Chem.MolFromMolFile(mol_in)
print('Adding hydrogens')
rdkitmolh = Chem.AddHs(rdkitmol, addCoords=True)
# ensure the chiral centers are all defined
Chem.AssignAtomChiralTagsFromStructure(rdkitmolh)
ligand_mol = Molecule(rdkitmolh)

print('Reading protein')
protein_pdb = PDBFile(pdb_in)

# Use Modeller to combine the protein and ligand into a complex
print('Preparing complex')
modeller = Modeller(protein_pdb.topology, protein_pdb.positions)

# This next bit is black magic.
# Modeller needs topology and positions. Lots of trial and error found that this is what works to get these from
# an openforcefield Molecule object that was created from a RDKit molecule.
# The topology part is described in the openforcefield API but the positions part grabs the first (and only)
# conformer and passes it to Modeller. It works. Don't ask why!
modeller.add(ligand_mol.to_topology().to_openmm(), ligand_mol.conformers[0])
# modeller.topology.setPeriodicBoxVectors(
#     [Vec3(x=8.461, y=0.0, z=0.0),
#     Vec3(x=0.0, y=8.461, z=0.0),
#     Vec3(x=0.0, y=0.0, z=8.461)])

print('System has %d atoms' % modeller.topology.getNumAtoms())

with open(output_complex, 'w') as outfile:
    PDBFile.writeFile(modeller.topology, modeller.positions, outfile)

# Initialize a SystemGenerator using the GAFF for the ligand
print('Preparing system')
forcefield_kwargs = {'constraints': app.HBonds, 'rigidWater': True, 'removeCMMotion': False, 'hydrogenMass': 4*unit.amu}
system_generator = SystemGenerator(
    forcefields=['amber/ff14SB.xml'],
    small_molecule_forcefield='gaff-2.11',
    forcefield_kwargs=forcefield_kwargs)

system = system_generator.create_system(modeller.topology, molecules=ligand_mol)

friction_coeff = 1 / unit.picosecond
step_size = 0.002 * unit.picoseconds
duration = (step_size * num_steps).value_in_unit(unit.nanoseconds)
print('Simulating for {} ns'.format(duration))

integrator = LangevinIntegrator(temperature, friction_coeff, step_size)
# system.addForce(openmm.MonteCarloBarostat(1*unit.atmospheres, temperature, 25))
print('Uses Periodic box: {}, Default Periodic box: {}'.format(
    system.usesPeriodicBoundaryConditions(), system.getDefaultPeriodicBoxVectors()))

simulation = Simulation(modeller.topology, system, integrator, platform=platform)
simulation.context.setPositions(modeller.positions)
print('Minimising ...')
simulation.minimizeEnergy()

# write out the minimised PDB
with open(output_min, 'w') as outfile:
    PDBFile.writeFile(modeller.topology, simulation.context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(), file=outfile, keepIds=True)

# equilibrate
simulation.context.setVelocitiesToTemperature(temperature)
print('Equilibrating ...')
simulation.step(equilibration_steps)

# Run the simulation.
# The enforcePeriodicBox arg to the reporters is important.
# It's a bit counter-intuitive that the value needs to be False, but this is needed to ensure that
# all parts of the simulation end up in the same periodic box when being output.
# simulation.reporters.append(PDBReporter(output_traj_pdb, reporting_interval, enforcePeriodicBox=False))
simulation.reporters.append(DCDReporter(output_traj_dcd, reporting_interval, enforcePeriodicBox=False))
simulation.reporters.append(StateDataReporter(sys.stdout, reporting_interval * 5, step=True, potentialEnergy=True, temperature=True))
print('Starting simulation with', num_steps, 'steps ...')
t1 = time.time()
simulation.step(num_steps)
t2 = time.time()
print('Simulation complete in {} mins at {}. Total wall clock time was {} mins'.format(
    round((t2 - t1) / 60, 3), temperature, round((t2 - t0) / 60, 3)))
print('Simulation time was', round(duration, 3), 'ns')
