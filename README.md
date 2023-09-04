# OpenMM simulation example

This repo is an effort to work out how to run a simulation of a protein-ligand complex with OpenMM.

What I'm looking for is for this to be as simple as possible:
1. read protein from PDB file
2. read ligand ideally from molfile or SDF
3. combine the protein and ligand into a complex
4. parameterise the ligand using GAFF
5. simulate
6. analyse

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
* better handling of commandline options (using argparse)

## simulateProtein.py: Basic Protein Simulation

Minimal example of setting up a protein and doing a simple minimisation.

Not much to say about this.

## simulateComplex.py: Basic Complex Simulation

How to set up a protein-ligand complex simulation.

After some help from @jchodera I put together this example to help illustrate this.

```
python simulateComplex.py -p protein.pdb -l ligand1.mol -o output -s 5000
```
Options:
```
$ python simulateComplex.py --help
usage: simulateComplex.py [-h] -p PROTEIN -l LIGAND [-o OUTPUT] [-s STEPS] [-z STEP_SIZE] [-f FRICTION_COEFF] [-i INTERVAL] [-t TEMPERATURE] [-e EQUILIBRATION_STEPS] [--protein-force-field PROTEIN_FORCE_FIELD]
                          [--ligand-force-field LIGAND_FORCE_FIELD]

simulateComplex

options:
  -h, --help            show this help message and exit
  -p PROTEIN, --protein PROTEIN
                        Protein PDB file
  -l LIGAND, --ligand LIGAND
                        Ligand molfile
  -o OUTPUT, --output OUTPUT
                        Base name for output files (output.dcd
  -s STEPS, --steps STEPS
                        Number of steps
  -z STEP_SIZE, --step-size STEP_SIZE
                        Step size (ps
  -f FRICTION_COEFF, --friction-coeff FRICTION_COEFF
                        Friction coefficient (ps
  -i INTERVAL, --interval INTERVAL
                        Reporting interval
  -t TEMPERATURE, --temperature TEMPERATURE
                        Temperature (K)
  -e EQUILIBRATION_STEPS, --equilibration-steps EQUILIBRATION_STEPS
                        Number of equilibration steps
  --protein-force-field PROTEIN_FORCE_FIELD
                        Protein force field
  --ligand-force-field LIGAND_FORCE_FIELD
                        Ligand force field
```

The protein is then read in PDB format and added to a Modeller object.
The ligand is then added to the Modeller to generate the complex.
The system is then prepared using the appropriate force fields.

Try the simulation as:

```
$ python simulateComplex.py -p protein.pdb -l ligand1.mol -o output -s 5000
simulateComplex:  Namespace(protein='protein.pdb', ligand='ligand1.mol', output='output', steps=5000, step_size=0.002, friction_coeff=1, interval=1000, temperature=300, equilibration_steps=200, protein_force_field='amber/ff14SB.xml', ligand_force_field='gaff-2.11')
Processing protein.pdb and ligand1.mol with 5000 steps generating outputs output_complex.pdb output_minimised.pdb output_traj.pdb output_traj.dcd
Using platform CUDA
Set precision for platform CUDA to mixed
Reading protein
Preparing complex
System has 4666 atoms
Preparing system
Simulating for 0.01 ns
Minimising ...
Equilibrating ...
Starting simulation with 5000 steps ...
#"Step","Potential Energy (kJ/mole)","Temperature (K)"
5000,-19712.275550323426,295.415879346069
Simulation complete in 0.02 mins at 300 K. Total wall clock time was 0.103 mins
Simulation time was 0.01 ns
```

The files `output_complex.pdb`, `output_minimised.pdb` and `output_traj.dcd` are generated.
The first is the complex, the second is that complex mimimised ready for simulation, the third the MD trajectory in DCD format.

See the code for details and gotchas.

## simulateComplexWithSolvent.py: Simulation with explicit solvent

```
python simulateComplexWithSolvent.py -p protein.pdb -l ligand1.mol -o output -s 5000
```

Build on the previous [simulateComplex.py]() example by including explicit solvent and a periodic box.
The system now has 58,052 atoms and takes quite a lot longer to simulate.
On my laptop's GeForce GTX 1050 GPU it takes almost 2 mins. A 1ns simulation takes just under one hour. (these were using
OpenMM 7.4).
On a newer Geforce RTX 3060 GPU (a modern mid-range gamers card) it takes just over half a minute, and the 1ns simulation
takes approx 15 mins.

Output is similar to the previous example.
See the code for details and gotchas.


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

The MD trajectories are analysed using [MDTraj](http://mdtraj.org/) and the script [analyse.py]().
```
python analyse.py trajectory.dcd topology.pdb output
```
This requires the trajectory to be written out using the DCD reporter. The topology can be read from the minimised
starting point of the MD run. This can be used for simulations with or without water.

The RMSD of the ligand and the protein C-alpha atoms compared to the start of the trajectory are displayed in a chart
that is generated using [Plotly](https://plotly.com/graphing-libraries/)) with the name output.svg (where 'output' is the
last parameter passed to the `analyse.py` script).

Example analysis:
![Example analysis](analyse.svg?raw=true "Example analysis]")

For complexes that are stable the RMSDs should not change dramatically. For a complex that is unstable the ligand may 
detach from the protein and the RMSD will increase dramatically. Relatively long simulations will be needed, maybe in the 
order of 100s of ns (input welcome on this and on how valid these simulations will be without explicit water).

## Improvements

Suggestions for how to improve these scripts and/or additional examples are welcome.
