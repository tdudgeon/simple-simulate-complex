# OpenMM simulation example

This repo is an effort to work out how to run a simulation of a protein-ligand complex with OpenMM.

What I'm looking for is for this to be as simple as possible:
1. read protein from PDB file
2. read ligand ideally from molfile or SDF
3. combine the protein and ligand into a complex
4. parameterise the ligand using GAFF
5. simulate
6. analyse

You might also want to look at this detailed example from openforcefields:
https://github.com/openforcefield/openff-toolkit/tree/stable/examples/toolkit_showcase.
It is a more elegant approach as it would allow to prepare a system for multiple MM toolkits (Gromacs, Amber etc.)
However, I can't get it to work for a number of reasons, so the examples here stick to using mostly OpenMM tooling.

Many thanks to @jchodera and others in the OpenMM community for help in putting these scripts together.

## Versions and branches

OpenMM and its related tools change quite a bit over time. The original work here was based on OpenMM version 7.4,
and the tools have now been updated for newer versions. Because the code is somewhat version dependent I have created
branches for specific versions. The code on master branch will typically be for the latest version I have used (not
always the latest version of OpenMM).

### History

#### 2020
* initial version based on OpenMM 7.4
* main tools are for simulating a ligand-protein complex with and without solvent

#### Sep 2023
* updated for OpenMM 7.7
* combined `simulateComplex.py` and `simulateComplexWithSolvent.py` into a single script
* better handling of commandline options (using argparse)
* many more options specifiable

## simulateProtein.py: Basic Protein Simulation

```
python simulateProtein.py protein.pdb
```

Minimal example of setting up a protein and doing a simple minimisation.

Not much to say about this.

## simulateComplex.py: Protein-ligand Complex Simulation

How to set up a protein-ligand complex simulation.
This script replaces the original one and `simulateComplexWithSolvent.py` and provides a single script that assembles
a protein and a ligand and optionally allows to solvate the system.

Options:
```
$ python simulateComplex.py -h
usage: simulateComplex.py [-h] -p PROTEIN -l LIGAND [-o OUTPUT] [-s STEPS] [-z STEP_SIZE] [-f FRICTION_COEFF] [-i INTERVAL] [-t TEMPERATURE] [--solvate] [--padding PADDING]
                          [--water-model {tip3p,spce,tip4pew,tip5p,swm4ndp}] [--positive-ion POSITIVE_ION] [--negative-ion NEGATIVE_ION] [--ionic-strength IONIC_STRENGTH] [--no-neutralize] [-e EQUILIBRATION_STEPS]
                          [--protein-force-field PROTEIN_FORCE_FIELD] [--ligand-force-field LIGAND_FORCE_FIELD] [--water-force-field WATER_FORCE_FIELD]

simulateComplexWithSolvent

options:
  -h, --help            show this help message and exit
  -p PROTEIN, --protein PROTEIN
                        Protein PDB file
  -l LIGAND, --ligand LIGAND
                        Ligand molfile
  -o OUTPUT, --output OUTPUT
                        Base name for output files
  -s STEPS, --steps STEPS
                        Number of steps
  -z STEP_SIZE, --step-size STEP_SIZE
                        Step size (ps
  -f FRICTION_COEFF, --friction-coeff FRICTION_COEFF
                        Friction coefficient (ps)
  -i INTERVAL, --interval INTERVAL
                        Reporting interval
  -t TEMPERATURE, --temperature TEMPERATURE
                        Temperature (K)
  --solvate             Add solvent box
  --padding PADDING     Padding for solvent box (A)
  --water-model {tip3p,spce,tip4pew,tip5p,swm4ndp}
                        Water model for solvation
  --positive-ion POSITIVE_ION
                        Positive ion for solvation
  --negative-ion NEGATIVE_ION
                        Negative ion for solvation
  --ionic-strength IONIC_STRENGTH
                        Ionic strength for solvation
  --no-neutralize       Don't add ions to neutralize
  -e EQUILIBRATION_STEPS, --equilibration-steps EQUILIBRATION_STEPS
                        Number of equilibration steps
  --protein-force-field PROTEIN_FORCE_FIELD
                        Protein force field
  --ligand-force-field LIGAND_FORCE_FIELD
                        Ligand force field
  --water-force-field WATER_FORCE_FIELD
                        Ligand force field
```

Many of the options are related to the solvation and are not needed unless you use the `--solvate` option. Most options
have sensible defaults.

The protein is read in PDB format and added to a Modeller object.
The ligand is then added to the Modeller to generate the complex.
If `--solvate` is specified a solvent box is added.
The system is then prepared using the appropriate force fields, the complex minimised, then equilibrated and finally
the MD simulation is run.

Try the simulation without solvation as:

```
$ python simulateComplex.py -p protein.pdb -l ligand1.mol
simulateComplexWithSolvent:  Namespace(protein='protein.pdb', ligand='ligand1.mol', output='output', steps=5000, step_size=0.002, friction_coeff=1, interval=1000, temperature=300, solvate=False, padding=10, water_model='tip3p', positive_ion='Na+', negative_ion='Cl-', ionic_strength=0.0, no_neutralize=False, equilibration_steps=200, protein_force_field='amber/ff14SB.xml', ligand_force_field='gaff-2.11', water_force_field='amber/tip3p_standard.xml')
Processing protein.pdb and ligand1.mol with 5000 steps generating outputs output_complex.pdb output_minimised.pdb output_traj.dcd
Using platform CUDA
Set precision for platform CUDA to mixed
Reading ligand
Preparing system
Reading protein
Preparing complex
System has 4645 atoms
Adding ligand...
System has 4666 atoms
Simulating for 0.01 ns
No Periodic Box
Minimising ...
Equilibrating ...
Starting simulation with 5000 steps ...
#"Step","Potential Energy (kJ/mole)","Temperature (K)"
5000,-19906.630997571585,297.9683027771655
Simulation complete in 0.02 mins at 300 K. Total wall clock time was 0.106 mins
Simulation time was 0.01 ns
```

The files `output_complex.pdb`, `output_minimised.pdb` and `output_traj.dcd` are generated.
The first is the complex, the second is that complex mimimised ready for simulation, the third the MD trajectory in DCD
format.

See the code for details and gotchas.

To run with solvation use something like this:

```
$ python simulateComplex.py -p protein.pdb -l ligand1.mol --solvate
simulateComplexWithSolvent:  Namespace(protein='protein.pdb', ligand='ligand1.mol', output='output', steps=5000, step_size=0.002, friction_coeff=1, interval=1000, temperature=300, solvate=True, padding=10, water_model='tip3p', positive_ion='Na+', negative_ion='Cl-', ionic_strength=0.0, no_neutralize=False, equilibration_steps=200, protein_force_field='amber/ff14SB.xml', ligand_force_field='gaff-2.11', water_force_field='amber/tip3p_standard.xml')
Processing protein.pdb and ligand1.mol with 5000 steps generating outputs output_complex.pdb output_minimised.pdb output_traj.dcd
Using platform CUDA
Set precision for platform CUDA to mixed
Reading ligand
Preparing system
Reading protein
Preparing complex
System has 4645 atoms
Adding ligand...
System has 4666 atoms
Adding solvent...
System has 58052 atoms
Simulating for 0.01 ns
Default Periodic box: [Quantity(value=Vec3(x=8.481300000000001, y=0.0, z=0.0), unit=nanometer), Quantity(value=Vec3(x=0.0, y=8.481300000000001, z=0.0), unit=nanometer), Quantity(value=Vec3(x=0.0, y=0.0, z=8.481300000000001), unit=nanometer)]
Minimising ...
Equilibrating ...
Starting simulation with 5000 steps ...
#"Step","Potential Energy (kJ/mole)","Temperature (K)"
5000,-742769.3383788001,301.45460839045137
Simulation complete in 0.145 mins at 300 K. Total wall clock time was 0.508 mins
Simulation time was 0.01 ns
```

The system now has 58,052 atoms and takes quite a lot longer to simulate.
On my laptop's GeForce GTX 1050 GPU it takes almost 2 mins. A 1ns simulation takes just under one hour. (these were using
OpenMM 7.4).
On a newer Geforce RTX 3060 GPU (a modern mid-range gamers card) it takes just over half a minute, and the 1ns simulation
takes approx 15 mins, a 100ns simulation about 20 hours.

Output is similar to the previous example.

Note that when adding solvent you are creating a periodic box system, and this can introduce weird visual quirks.
See [this discussion](https://github.com/openmm/openmm/issues/4218) for more details.
Those quirks can be resolved by using the `analyse.py` script that is described below.

## Protein and ligand preparation

*Note: these scripts have not yet been updated for OpenMM 77*

The previous methods were a bit of a cheat as they used a protein that had been fully prepared for
simulation. It's pretty unlikely you will start with a protein in that state. There are 2 scripts that
illustrate how preparation can be done. The aim is to be able to do this entirely within OpenMM, but it seems
that's not quite possible.

The scripts are:

[prepareProtein.py]()
```
 python prepareProtein.py protein.pdb protein
```
This strips out everything that is not the protein, fixes problems in the protein, adds hydrogens and writes the
file `protein_prepared.pdb`. That file can be used as inputs to the previous simulations.

[prepareComplex.py]()
```
 python prepareComplex.py data/Mpro-x0387_0_apo-desolv.pdb Mpro-x0387_0
```
This aims to build a PDB file with protein and ligand (and optionally the crystallographic waters) that is
ready for simulation. It writes the files `protein_prepared.pdb` and `ligand_prepared.pdb`.
It doesn't do everything that's needed, so other toolkits will be required:
- ligand does not have hydrogens added
- ligand can only be written to PDB format

## Analysis

The MD trajectories are analysed using the script [analyse.py]() which uses [MDTraj](http://mdtraj.org/).

Usage:
```
$ python analyse.py -h
usage: analyse.py [-h] -p PROTEIN -t TRAJECTORY -o OUTPUT [-r]

analyse

options:
  -h, --help            show this help message and exit
  -p PROTEIN, --protein PROTEIN
                        Protein PDB file
  -t TRAJECTORY, --trajectory TRAJECTORY
                        Trajectory DCD file
  -o OUTPUT, --output OUTPUT
                        Output base name
  -r, --remove-waters   Remove waters, salts etc.
```
This requires the trajectory to be written out using the DCD reporter. The topology can be read from the minimised
starting point of the MD run. This can be used for simulations with or without water.

The RMSD of the ligand and the protein C-alpha atoms compared to the start of the trajectory are displayed in a chart
that is generated using [Plotly](https://plotly.com/graphing-libraries/)) with the name output.svg.

The trajectory is also re-imaged so that the ligand is in the correct location with respect to the protein (the periodic
box can cause some wierd visual aberations), and the waters, ions etc. removed if required.

Example:
```
$ python analyse.py -p output_minimised.pdb -t output_traj.dcd -o output_reimaged -r
analyse:  Namespace(protein='output_minimised.pdb', trajectory='output_traj.dcd', output='output_reimaged', remove_waters=True)
Reading trajectory output_traj.dcd
Removing waters
Realigning
Writing re-imaged PDB output_reimaged.pdb
Writing re-imaged trajectory output_reimaged.dcd
Number of frames: 500
21 ligand atoms
1216 backbone atoms
Writing RMSD output to output_reimaged.svg
```

In this case 3 files are created:
* output_reimaged.svg - SVG of the RMSD of the backbone and ligand
* output_reimaged.pdb - the re-imaged PDB file
* output_reimaged.dcd - the re-imaged trajectory

Example RMSD analysis:
![Example analysis](analyse.svg?raw=true "Example analysis]")

The trajectory can be nicely viewed in [NGLView](http://nglviewer.org/ngl/).
First load the PDB file, then click on the menu for the item in the selector on the right, find the Trajectory section
and load the DCD file.

For complexes that are stable the RMSDs should not change dramatically. For a complex that is unstable the ligand may 
detach from the protein and the RMSD will increase dramatically. Relatively long simulations will be needed, maybe in the 
order of 100s of ns (input welcome on this and on how valid these simulations will be without explicit water).

## Improvements

Suggestions for how to improve these scripts and/or additional examples are welcome.
