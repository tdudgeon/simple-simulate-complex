import sys, time
from openforcefield.topology import Molecule
from openmmforcefields.generators import SystemGenerator
from simtk import unit, openmm
from simtk.openmm import app, Platform, LangevinIntegrator
from simtk.openmm.app import PDBFile, Simulation, Modeller, PDBReporter, StateDataReporter, DCDReporter
from rdkit import Chem


temperature = 300 * unit.kelvin
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


# Read the molfile into RDKit, add Hs and create an openforcefield Molecule object
print('Reading ligand')
rdkitmol = Chem.MolFromMolFile(mol_in)
print('Adding hydrogens')
rdkitmolh = Chem.AddHs(rdkitmol, addCoords=True)
# ensure the chiral centers are all defined
Chem.AssignAtomChiralTagsFromStructure(rdkitmolh)
ligand_mol = Molecule(rdkitmolh)

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
protein_pdb = PDBFile(pdb_in)

print('Preparing complex')
modeller = Modeller(protein_pdb.topology, protein_pdb.positions)
print('System has %d atoms' % modeller.topology.getNumAtoms())

# This next bit is black magic.
# Modeller needs topology and positions. Lots of trial and error found that this is what works to get these from
# an openforcefield Molecule object that was created from a RDKit molecule.
# The topology part is described in the openforcefield API but the positions part grabs the first (and only)
# conformer and passes it to Modeller. It works. Don't ask why!
modeller.add(ligand_mol.to_topology().to_openmm(), ligand_mol.conformers[0])

print('System has %d atoms' % modeller.topology.getNumAtoms())

# Solvate
print('Adding solvent...')
# we use the 'padding' option to define the periodic box. The PDB file does not contain any
# unit cell information so we just create a box that has a 10A padding around the complex.
modeller.addSolvent(system_generator.forcefield, model='tip3p', padding=10.0*unit.angstroms)
print('System has %d atoms' % modeller.topology.getNumAtoms())

with open(output_complex, 'w') as outfile:
    PDBFile.writeFile(modeller.topology, modeller.positions, outfile)

# Create the system using the SystemGenerator
system = system_generator.create_system(modeller.topology, molecules=ligand_mol)
integrator = LangevinIntegrator(temperature, 1 / unit.picosecond, 0.002 * unit.picoseconds)
system.addForce(openmm.MonteCarloBarostat(1*unit.atmospheres, temperature, 25))
print('Default Periodic box:', system.getDefaultPeriodicBoxVectors())

simulation = Simulation(modeller.topology, system, integrator, platform=platform)
context = simulation.context
context.setPositions(modeller.positions)

print('Minimising ...')
simulation.minimizeEnergy()

# Write out the minimised PDB. The 'enforcePeriodicBox=False' bit is important otherwise the different
# components can end up in different periodic boxes resulting in really strange looking output.
with open(output_min, 'w') as outfile:
    PDBFile.writeFile(modeller.topology, context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(), file=outfile, keepIds=True)

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
t0 = time.time()
simulation.step(num_steps)
t1 = time.time()
print('Simulation complete in', t1 - t0, 'seconds at', temperature)
