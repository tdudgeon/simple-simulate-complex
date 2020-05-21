# OpenMM simulation example

This repo is an effort to work out how to run a simulation of a protein-ligand complex with OpenMM.

What I'm looking for is for this to be as simple as possible:
1. read protein from PDB file
2. read ligand ideally from molfile, but can cope with PDB
3. combine the protein and ligand into a complex
4. parameterise the ligand using GAFF
5. simulate

## simulateComplex.py: Basic Simulation

How to set up a protein-ligand complex simulation.

After some help from @jchodera I put together this example
to help illustrate this. Note the issues that seem present at the bottom of this doc.

Note: currently this is running with OpenMM 7.4.1 because https://github.com/openmm/openmm/issues/2683
makes it difficult to use 3.5.0.

The example:

1. protein.pdb - a protein with hydrogens ready to simulate (as done in [](simulateProtein.py))
2. ligand1.mol/sdf/pdb - ligand with hydrogens in mol, sdf and pdb formats

Try the simulation as:

```
$ python simulateComplex.py 
Warning: Unable to load toolkit 'OpenEye Toolkit'. The Open Force Field Toolkit does not require the OpenEye Toolkits, and can use RDKit/AmberTools instead. However, if you have a valid license for the OpenEye Toolkits, consider installing them for faster performance and additional file format support: https://docs.eyesopen.com/toolkits/python/quickstart-python/linuxosx.html OpenEye offers free Toolkit licenses for academics: https://www.eyesopen.com/academic-licensing
Molecule with name 'CCOC(=O)c1ccccc1' and SMILES '[H][C]1=[C]([H])[C]([H])=[C]([C](=[O])[O][C]([H])([H])[C]([H])([H])[H])[C]([H])=[C]1[H]'
Reading protein
Reading ligand
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM    2  C   UNL     1       9.638  -2.277  28.027  1.00  0.00           C  , HETATM    1  C   UNL     1       9.204  -3.694  28.363  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM    4  C   UNL     1      10.656  -1.091  26.262  1.00  0.00           C  , HETATM    2  C   UNL     1       9.638  -2.277  28.027  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM    5  O   UNL     1      10.704  -0.055  26.905  1.00  0.00           O  , HETATM    3  O   UNL     1      10.180  -2.283  26.700  1.00  0.00           O  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM    6  C   UNL     1      11.242  -1.246  24.906  1.00  0.00           C  , HETATM    4  C   UNL     1      10.656  -1.091  26.262  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM    7  C   UNL     1      11.968  -0.168  24.378  1.00  0.00           C  , HETATM    6  C   UNL     1      11.242  -1.246  24.906  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM    8  C   UNL     1      12.609  -0.292  23.145  1.00  0.00           C  , HETATM    7  C   UNL     1      11.968  -0.168  24.378  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM    9  C   UNL     1      12.525  -1.487  22.431  1.00  0.00           C  , HETATM    8  C   UNL     1      12.609  -0.292  23.145  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM   10  C   UNL     1      11.796  -2.560  22.943  1.00  0.00           C  , HETATM    9  C   UNL     1      12.525  -1.487  22.431  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM   11  C   UNL     1      11.153  -2.441  24.177  1.00  0.00           C  , HETATM   10  C   UNL     1      11.796  -2.560  22.943  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM   13  H   UNL     1      10.050  -4.348  28.312  1.00  0.00           H  , HETATM   12  H   UNL     1       8.794  -3.717  29.351  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM   14  H   UNL     1       8.462  -4.016  27.662  1.00  0.00           H  , HETATM   13  H   UNL     1      10.050  -4.348  28.312  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM   15  H   UNL     1      10.383  -1.951  28.722  1.00  0.00           H  , HETATM   14  H   UNL     1       8.462  -4.016  27.662  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM   16  H   UNL     1       8.804  -1.609  28.085  1.00  0.00           H  , HETATM   15  H   UNL     1      10.383  -1.951  28.722  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM   17  H   UNL     1      12.027   0.719  24.902  1.00  0.00           H  , HETATM   16  H   UNL     1       8.804  -1.609  28.085  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM   18  H   UNL     1      13.146   0.502  22.761  1.00  0.00           H  , HETATM   17  H   UNL     1      12.027   0.719  24.902  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM   19  H   UNL     1      13.004  -1.577  21.521  1.00  0.00           H  , HETATM   18  H   UNL     1      13.146   0.502  22.761  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM   20  H   UNL     1      11.732  -3.442  22.411  1.00  0.00           H  , HETATM   19  H   UNL     1      13.004  -1.577  21.521  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM   21  H   UNL     1      10.610  -3.235  24.553  1.00  0.00           H  , HETATM   20  H   UNL     1      11.732  -3.442  22.411  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
Preparing complex
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4648  C   UNL B   1       9.638  -2.277  28.027  1.00  0.00           C  , HETATM 4647  C   UNL B   1       9.204  -3.694  28.363  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4650  C   UNL B   1      10.656  -1.091  26.262  1.00  0.00           C  , HETATM 4648  C   UNL B   1       9.638  -2.277  28.027  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4651  O   UNL B   1      10.704  -0.055  26.905  1.00  0.00           O  , HETATM 4649  O   UNL B   1      10.180  -2.283  26.700  1.00  0.00           O  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4652  C   UNL B   1      11.242  -1.246  24.906  1.00  0.00           C  , HETATM 4650  C   UNL B   1      10.656  -1.091  26.262  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4653  C   UNL B   1      11.968  -0.168  24.378  1.00  0.00           C  , HETATM 4652  C   UNL B   1      11.242  -1.246  24.906  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4654  C   UNL B   1      12.609  -0.292  23.145  1.00  0.00           C  , HETATM 4653  C   UNL B   1      11.968  -0.168  24.378  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4655  C   UNL B   1      12.525  -1.487  22.431  1.00  0.00           C  , HETATM 4654  C   UNL B   1      12.609  -0.292  23.145  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4656  C   UNL B   1      11.796  -2.560  22.943  1.00  0.00           C  , HETATM 4655  C   UNL B   1      12.525  -1.487  22.431  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4657  C   UNL B   1      11.153  -2.441  24.177  1.00  0.00           C  , HETATM 4656  C   UNL B   1      11.796  -2.560  22.943  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4659  H   UNL B   1      10.050  -4.348  28.312  1.00  0.00           H  , HETATM 4658  H   UNL B   1       8.794  -3.717  29.351  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4660  H   UNL B   1       8.462  -4.016  27.662  1.00  0.00           H  , HETATM 4659  H   UNL B   1      10.050  -4.348  28.312  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4661  H   UNL B   1      10.383  -1.951  28.722  1.00  0.00           H  , HETATM 4660  H   UNL B   1       8.462  -4.016  27.662  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4662  H   UNL B   1       8.804  -1.609  28.085  1.00  0.00           H  , HETATM 4661  H   UNL B   1      10.383  -1.951  28.722  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4663  H   UNL B   1      12.027   0.719  24.902  1.00  0.00           H  , HETATM 4662  H   UNL B   1       8.804  -1.609  28.085  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4664  H   UNL B   1      13.146   0.502  22.761  1.00  0.00           H  , HETATM 4663  H   UNL B   1      12.027   0.719  24.902  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4665  H   UNL B   1      13.004  -1.577  21.521  1.00  0.00           H  , HETATM 4664  H   UNL B   1      13.146   0.502  22.761  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4666  H   UNL B   1      11.732  -3.442  22.411  1.00  0.00           H  , HETATM 4665  H   UNL B   1      13.004  -1.577  21.521  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4667  H   UNL B   1      10.610  -3.235  24.553  1.00  0.00           H  , HETATM 4666  H   UNL B   1      11.732  -3.442  22.411  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
Preparing system
Warning: In AmberToolsToolkitwrapper.compute_partial_charges_am1bcc: Molecule 'CCOC(=O)c1ccccc1' has more than one conformer, but this function will only generate charges for the first one.

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

Info: Total number of electrons: 80; net charge: 0

Running: /home/timbo/miniconda3/envs/openmm/bin/sqm -O -i sqm.in -o sqm.out


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


Minimising
Starting simulation
#"Step","Potential Energy (kJ/mole)","Temperature (K)"
1000,-21518.453125,267.3147547214519
2000,-19785.859375,298.57485160064016
3000,-19828.603515625,306.5371419137492
4000,-19524.287109375,298.33132752385455
5000,-19886.439453125,298.03124036938266
6000,-20094.91015625,302.5826579116727
7000,-19967.19921875,299.14510478741903
8000,-20275.4765625,300.6640170468736
9000,-20103.796875,295.850503453827
10000,-20001.4921875,295.04521624908136
Done
```

The file `output.pdb` is generated and does contain a trajectory of the complex.
So that example is now working, but I don't see it as optimal because of these issues:

1. The ligand has to be read twice, once as SDF so that Amber can parameterise it, once as
PDB so that if can be combined into the complex. It should be possible to read the molecule
just once in any common format (e.g. molfile, SDF, PDB, MOL2, RDKit RWMol object ...).
2. When reading the ligand PDB a shed load of `WARNING: duplicate atom` warnings are issued.
And again when reading the complex.
Clearly these warnings are incorrect. The atoms are clearly not duplicates.
3. Amber seems to issue this warning
`Warning: In AmberToolsToolkitwrapper.compute_partial_charges_am1bcc: Molecule 'CCOC(=O)c1ccccc1' has more than one conformer, but this function will only generate charges for the first one.`
when reading the ligand. Clearly it is wrong. There is only one conformer present.
Probably nothing OpenMM can do about this.


## simulateComplexWithSolvent.py: Simulation with explicit solvent

```
python simulateComplexWithSolvent.py
```

Build on the previous [simulateComplex.py]() example by including explicit solvent.
This is not currently working because the minimisation step blows up (the ligand and protein move outside
the box of waters). Probably this is something related to the periodic box setup.

Another issue is defining the size of the box. This is currently hard-coded to values that look about 
right, but it needs to be calculated on the fly. Presumably OpenMM provides an easy way to do this,
but I can't find it!
