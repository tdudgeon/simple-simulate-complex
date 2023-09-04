import sys, time, argparse

from openff.toolkit import Molecule
from openmmforcefields.generators import SystemGenerator
import openmm
from openmm import app, unit, LangevinIntegrator
from openmm.app import PDBFile, Simulation, Modeller, PDBReporter, StateDataReporter, DCDReporter

import utils

t0 = time.time()


parser = argparse.ArgumentParser(description="simulateComplexWithSolvent")

parser.add_argument("-p", "--protein", required=True, help="Protein PDB file")
parser.add_argument("-l", "--ligand", required=True, help="Ligand molfile")
parser.add_argument("-o", "--output", default='output', help="Base name for output files")
parser.add_argument("-s", "--steps", type=int, default=5000, help="Number of steps")
parser.add_argument("-z", "--step-size", type=float, default=0.002, help="Step size (ps")
parser.add_argument("-f", "--friction-coeff", type=float, default=1, help="Friction coefficient (ps)")
parser.add_argument("-i", "--interval", type=int, default=1000, help="Reporting interval")
parser.add_argument("-t", "--temperature", type=int, default=300, help="Temperature (K)")
parser.add_argument("--padding", type=float, default=10, help="Padding for solvent box (A)")
parser.add_argument("-e", "--equilibration-steps", type=int, default=200, help="Number of equilibration steps")
parser.add_argument("--protein-force-field", default='amber/ff14SB.xml', help="Protein force field")
parser.add_argument("--ligand-force-field", default='gaff-2.11', help="Ligand force field")
parser.add_argument("--water-force-field", default='amber/tip3p_standard.xml', help="Ligand force field")

args = parser.parse_args()
print("simulateComplexWithSolvent: ", args)

pdb_in = args.protein
mol_in = args.ligand
output_base = args.output
output_complex = output_base + '_complex.pdb'
output_traj_pdb = output_base + '_traj.pdb'
output_traj_dcd = output_base + '_traj.dcd'
output_min = output_base + '_minimised.pdb'
num_steps = args.steps
reporting_interval = args.interval
temperature = args.temperature * unit.kelvin
equilibration_steps = args.equilibration_steps
print('Processing', pdb_in, 'and', mol_in, 'with', num_steps, 'steps generating outputs',
      output_complex, output_min, output_traj_pdb, output_traj_dcd)

# get the chosen or fastest platform
platform = utils.get_platform()

print('Reading ligand')
ligand_mol = Molecule.from_file(mol_in)

print('Preparing system')
# Initialize a SystemGenerator using the GAFF for the ligand and tip3p for the water.
forcefield_kwargs = {'constraints': app.HBonds, 'rigidWater': True, 'removeCMMotion': False, 'hydrogenMass': 4*unit.amu }
system_generator = SystemGenerator(
    forcefields=[args.protein_force_field, args.water_force_field],
    small_molecule_forcefield=args.ligand_force_field,
    molecules=[ligand_mol],
    forcefield_kwargs=forcefield_kwargs)

# Use Modeller to combine the protein and ligand into a complex
print('Reading protein')
protein_pdb = PDBFile(pdb_in)

print('Preparing complex')
# The topology is described in the openforcefield API

modeller = Modeller(protein_pdb.topology, protein_pdb.positions)
print('System has %d atoms' % modeller.topology.getNumAtoms())

# The topology is described in the openforcefield API
print('Adding ligand...')
lig_top = ligand_mol.to_topology()
modeller.add(lig_top.to_openmm(), lig_top.get_positions().to_openmm())

print('System has %d atoms' % modeller.topology.getNumAtoms())

# Solvate
print('Adding solvent...')
# we use the 'padding' option to define the periodic box. The PDB file does not contain any
# unit cell information so we just create a box that has a 10A (default) padding around the complex.
modeller.addSolvent(system_generator.forcefield, model='tip3p', padding=args.padding * unit.angstroms)
print('System has %d atoms' % modeller.topology.getNumAtoms())

with open(output_complex, 'w') as outfile:
    PDBFile.writeFile(modeller.topology, modeller.positions, outfile)

# Create the system using the SystemGenerator
system = system_generator.create_system(modeller.topology, molecules=ligand_mol)

friction_coeff = args.friction_coeff / unit.picosecond
step_size = args.step_size * unit.picoseconds
duration = (step_size * num_steps).value_in_unit(unit.nanoseconds)
print('Simulating for {} ns'.format(duration))

integrator = LangevinIntegrator(temperature, friction_coeff, step_size)
system.addForce(openmm.MonteCarloBarostat(1 * unit.atmospheres, temperature, 25))
print('Uses Periodic box: {}, Default Periodic box: {}'.format(
    system.usesPeriodicBoundaryConditions(), system.getDefaultPeriodicBoxVectors()))

simulation = Simulation(modeller.topology, system, integrator, platform=platform)
context = simulation.context
context.setPositions(modeller.positions)

print('Minimising ...')
simulation.minimizeEnergy()

# Write out the minimised PDB.
# Note: the enforcePeriodicBox=True option does not seem to work correctly
with open(output_min, 'w') as outfile:
    PDBFile.writeFile(modeller.topology, context.getState(getPositions=True, enforcePeriodicBox=True).getPositions(), file=outfile, keepIds=True)

# equilibrate
simulation.context.setVelocitiesToTemperature(temperature)
print('Equilibrating ...')
simulation.step(equilibration_steps)

# Run the simulation.
# The enforcePeriodicBox=True options does not seem to work correctly.
# simulation.reporters.append(PDBReporter(output_traj_pdb, reporting_interval, enforcePeriodicBox=True))
simulation.reporters.append(DCDReporter(output_traj_dcd, reporting_interval, enforcePeriodicBox=True))
simulation.reporters.append(StateDataReporter(sys.stdout, reporting_interval * 5, step=True, potentialEnergy=True, temperature=True))
print('Starting simulation with', num_steps, 'steps ...')
t1 = time.time()
simulation.step(num_steps)
t2 = time.time()
print('Simulation complete in {} mins at {}. Total wall clock time was {} mins'.format(
    round((t2 - t1) / 60, 3), temperature, round((t2 - t0) / 60, 3)))
print('Simulation time was', round(duration, 3), 'ns')
