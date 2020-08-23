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

## simulateProtein.py: Basic Protein Simulation

Minimal example of setting up a protein and doing a simple minimisation.

Not much to say about this.

## simulateComplex.py: Basic Complex Simulation

How to set up a protein-ligand complex simulation.

After some help from @jchodera I put together this example to help illustrate this.

Note: currently this is running with OpenMM 7.4.2 because https://github.com/openmm/openmm/issues/2683
makes it difficult to use 7.5.0.

```
 python simulateComplex.py protein.pdb ligand.mol output 5000
```
The arguments:

1. protein.pdb - a protein with hydrogens ready to simulate (as done in [simulateProtein.py]()
2. ligand.mol - ligand in molfile format
3. The name to use as the base name of the results.
4. The number of iterations for the simulation.

The ligand is read using RDKit and then processed to:
* Add hydrogens
* Define the stereochemistry

The protein is then read in PDB format and added to a Modeller object.
The ligand is then added to the Modeller to generate the complex.
The system is then prepared using the appropriate force fields.


Try the simulation as:

```
$ python simulateComplex.py Mpro-x0387_0_fixed.pdb data/Mpro-x0387_0.mol Mpro-x0387_0 5000
Warning: Unable to load toolkit 'OpenEye Toolkit'. The Open Force Field Toolkit does not require the OpenEye Toolkits, and can use RDKit/AmberTools instead. However, if you have a valid license for the OpenEye Toolkits, consider installing them for faster performance and additional file format support: https://docs.eyesopen.com/toolkits/python/quickstart-python/linuxosx.html OpenEye offers free Toolkit licenses for academics: https://www.eyesopen.com/academic-licensing
Processing Mpro-x0387_0_fixed.pdb and data/Mpro-x0387_0.mol with 5000 steps generating outputs Mpro-x0387_0_complex.pdb Mpro-x0387_0_minimised.pdb Mpro-x0387_0_traj.pdb Mpro-x0387_0_traj.dcd
Set precision for platform CUDA to mixed
Reading ligand
Adding hydrogens
Reading protein
Preparing complex
System has 4673 atoms
Preparing system
Warning: In AmberToolsToolkitwrapper.compute_partial_charges_am1bcc: Molecule '' has more than one conformer, but this function will only generate charges for the first one.

Welcome to antechamber 17.3: molecular input file processor.

acdoctor mode is on: check and diagnosis problems in the input file.
-- Check Format for sdf File --
   Status: pass
-- Check Unusual Elements --
   Status: pass
-- Check Open Valences --
   Status: pass
-- Check Geometry --
      for those bonded   
      for those not bonded   
   Status: pass
-- Check Weird Bonds --
   Status: pass
-- Check Number of Units --
   Status: pass
acdoctor mode has completed checking the input file.

Info: Total number of electrons: 106; net charge: 0

Running: /home/timbo/miniconda3/envs/openmm-74/bin/sqm -O -i sqm.in -o sqm.out


Welcome to antechamber 17.3: molecular input file processor.

acdoctor mode is on: check and diagnosis problems in the input file.
-- Check Format for mol2 File --
   Status: pass
-- Check Unusual Elements --
   Status: pass
-- Check Open Valences --
   Status: pass
-- Check Geometry --
      for those bonded   
      for those not bonded   
   Status: pass
-- Check Weird Bonds --
   Status: pass
-- Check Number of Units --
   Status: pass
acdoctor mode has completed checking the input file.


Default Periodic box: [Quantity(value=Vec3(x=11.215900000000001, y=0.0, z=0.0), unit=nanometer), Quantity(value=Vec3(x=0.0, y=5.2655, z=0.0), unit=nanometer), Quantity(value=Vec3(x=0.0, y=0.0, z=4.3365), unit=nanometer)]
Minimising ...
Equilibrating ...
Starting simulation with 5000 steps ...
#"Step","Potential Energy (kJ/mole)","Temperature (K)"
1000,-12905.476935093218,279.598239742249
2000,-12666.84060298634,296.43099174154327
3000,-13110.495784052066,304.59169315932655
4000,-13014.61967772014,307.1202441925517
5000,-13061.841891369622,294.2178513158892
Simulation complete in 6.061023712158203 seconds at 300 K
```

The files `Mpro-x0387_0_complex.pdb`, `Mpro-x0387_0_minimised.pdb` and `Mpro-x0387_0_traj.dcd` are generated.
The first is the complex, the second is that complex mimimised ready for simulation, the third the MD trajectory in DCD format.

See the code for details and gotchas.

## simulateComplexWithSolvent.py: Simulation with explicit solvent

```
python simulateComplexWithSolvent.py Mpro-x0387_0_fixed.pdb data/Mpro-x0387_0.mol Mpro-x0387_0 5000
```

Build on the previous [simulateComplex.py]() example by including explicit solvent.
The system now has 58,052 atoms and takes quite a lot longer to simulate, almost 2 mins using
my laptop's GeForce GTX 1050 GPU. A 1ns simulation takes just under one hour.

Output is similar to the previous example.
See the code for details and gotchas.


## Protein and ligand preparation

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