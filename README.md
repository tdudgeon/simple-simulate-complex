# OpenMM simulation example

This repo is an effort to work out how to run a simulation of a protein-ligand complex with OpenMM.

There is information on this scattered around the internet, but I can't get anything to work.
Seems like as well as OpenMM using openmmforcefields and parmed can play a part in this, but 
despite trying lots of permutations I couldn't get anything to work so I put together this example
to help illustrate the problem.

The example:

1. protein.pdb - a protein with hydrogens ready to simulate (as done in [](simulateProtein.py))
2. ligand1.mol/pdb - ligand with hydrogens in mol and pdb formats

Try the simulation as:
```
$ python simulateComplex.py 
Warning: Unable to load toolkit 'OpenEye Toolkit'. The Open Force Field Toolkit does not require the OpenEye Toolkits, and can use RDKit/AmberTools instead. However, if you have a valid license for the OpenEye Toolkits, consider installing them for faster performance and additional file format support: https://docs.eyesopen.com/toolkits/python/quickstart-python/linuxosx.html OpenEye offers free Toolkit licenses for academics: https://www.eyesopen.com/academic-licensing
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
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4648  C   UNL   305       9.638  -2.277  28.027  1.00  0.00           C  , HETATM 4647  C   UNL   305       9.204  -3.694  28.363  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4650  C   UNL   305      10.656  -1.091  26.262  1.00  0.00           C  , HETATM 4648  C   UNL   305       9.638  -2.277  28.027  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4651  O   UNL   305      10.704  -0.055  26.905  1.00  0.00           O  , HETATM 4649  O   UNL   305      10.180  -2.283  26.700  1.00  0.00           O  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4652  C   UNL   305      11.242  -1.246  24.906  1.00  0.00           C  , HETATM 4650  C   UNL   305      10.656  -1.091  26.262  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4653  C   UNL   305      11.968  -0.168  24.378  1.00  0.00           C  , HETATM 4652  C   UNL   305      11.242  -1.246  24.906  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4654  C   UNL   305      12.609  -0.292  23.145  1.00  0.00           C  , HETATM 4653  C   UNL   305      11.968  -0.168  24.378  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4655  C   UNL   305      12.525  -1.487  22.431  1.00  0.00           C  , HETATM 4654  C   UNL   305      12.609  -0.292  23.145  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4656  C   UNL   305      11.796  -2.560  22.943  1.00  0.00           C  , HETATM 4655  C   UNL   305      12.525  -1.487  22.431  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4657  C   UNL   305      11.153  -2.441  24.177  1.00  0.00           C  , HETATM 4656  C   UNL   305      11.796  -2.560  22.943  1.00  0.00           C  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4659  H   UNL   305      10.050  -4.348  28.312  1.00  0.00           H  , HETATM 4658  H   UNL   305       8.794  -3.717  29.351  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4660  H   UNL   305       8.462  -4.016  27.662  1.00  0.00           H  , HETATM 4659  H   UNL   305      10.050  -4.348  28.312  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4661  H   UNL   305      10.383  -1.951  28.722  1.00  0.00           H  , HETATM 4660  H   UNL   305       8.462  -4.016  27.662  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4662  H   UNL   305       8.804  -1.609  28.085  1.00  0.00           H  , HETATM 4661  H   UNL   305      10.383  -1.951  28.722  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4663  H   UNL   305      12.027   0.719  24.902  1.00  0.00           H  , HETATM 4662  H   UNL   305       8.804  -1.609  28.085  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4664  H   UNL   305      13.146   0.502  22.761  1.00  0.00           H  , HETATM 4663  H   UNL   305      12.027   0.719  24.902  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4665  H   UNL   305      13.004  -1.577  21.521  1.00  0.00           H  , HETATM 4664  H   UNL   305      13.146   0.502  22.761  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4666  H   UNL   305      11.732  -3.442  22.411  1.00  0.00           H  , HETATM 4665  H   UNL   305      13.004  -1.577  21.521  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/internal/pdbstructure.py:537: UserWarning: WARNING: duplicate atom (HETATM 4667  H   UNL   305      10.610  -3.235  24.553  1.00  0.00           H  , HETATM 4666  H   UNL   305      11.732  -3.442  22.411  1.00  0.00           H  )
  warnings.warn("WARNING: duplicate atom (%s, %s)" % (atom, old_atom._pdb_string(old_atom.serial_number, atom.alternate_location_indicator)))
Preparing system
Did not recognize residue UNL; did you forget to call .add_molecules() to add it?
Traceback (most recent call last):
  File "simulateComplex.py", line 42, in <module>
    system = system_generator.create_system(complex_pdb.topology, molecules=ligand_mol)
  File "/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/openmmforcefields/generators/system_generators.py", line 307, in create_system
    system = self.forcefield.createSystem(topology, **forcefield_kwargs)
  File "/home/timbo/miniconda3/envs/openmm/lib/python3.7/site-packages/simtk/openmm/app/forcefield.py", line 1145, in createSystem
    raise ValueError('No template found for residue %d (%s).  %s' % (res.index+1, res.name, _findMatchErrors(self, res)))
ValueError: No template found for residue 305 (UNL).  The set of atoms is similar to NTYR, but it is missing 2 atoms.
```

I've tried lots of different variations but haven't got any to work yet!

What I'm looking for is for this to be as simple as possible:
1. read protein from PDB file
2. read ligand ideally from molfile, but can cope with PDB
3. combine the protein and ligand into a complex
4. parameterise the ligand using GAFF
5. simulate