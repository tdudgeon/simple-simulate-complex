from simtk.openmm.app import *
from simtk.openmm import *
from pdbfixer import PDBFixer

if len(sys.argv) != 4:
    print('Usage: python prepareComplex.py input.pdb ligand.mol system')
    exit(1)

pdb_in = sys.argv[1]
ligand_in = sys.argv[2]
outname = sys.argv[3]

# This PDB file contains:
# - the protein (single chain)
# - a ligand
# - 3 DMSO molecules
# - A number of waters
# No hydrogens are present.
# The C-teminal THR residue is missing an oxygen atom.
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

# The following removes the DMS components and retains the ligand and waters.
# If instead we want to remove the ligand it will be easier to use:
# fixer.removeHeterogens(True) or fixer.removeHeterogens(False)
# True keeps the waters, False removes them leaving only the protein.
# See prepareProtein.py for this in action.
modeller = Modeller(fixer.topology, fixer.positions)
toDelete = []
for res in modeller.topology.residues():
    if res.name == 'DMS':
        toDelete.append(res)
        print('Deleting', res)
modeller.delete(toDelete)


with open(outname + '_receptor.pdb', 'w') as outfile:
    PDBFile.writeFile(modeller.topology, modeller.positions, file=outfile, keepIds=True)


# Now the ligand. We write it to a file in PDB format. It would be good to write to Molfile format but seems that
# OpenMM does not support that.
# Note: this may not be the best way to do this. Other toolkits might be better.
# Note: your ligand may not be named 'LIG'
# Note: modeller.addHydrogens() does not work for ligands. We'll need to use another toolkit such as OpenBabel to do this.
modeller = Modeller(fixer.topology, fixer.positions)
toDelete = []
for res in modeller.topology.residues():
    if res.name != 'LIG':
        toDelete.append(res)
modeller.delete(toDelete)

with open(outname + '_ligand.pdb', 'w') as outfile:
    PDBFile.writeFile(modeller.topology, modeller.positions, file=outfile, keepIds=True)

print('Done')