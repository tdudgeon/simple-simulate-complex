import sys, argparse
import mdtraj as md
import plotly.graph_objects as go


# Typical usage:
# python analyse.py -p output_minimised.pdb -t output_traj.dcd -o output_reimaged -r

parser = argparse.ArgumentParser(description="analyse")

parser.add_argument("-p", "--protein", required=True, help="Protein PDB file")
parser.add_argument("-t", "--trajectory", required=True, help="Trajectory DCD file")
parser.add_argument("-o", "--output", required=True, help="Output base name")
parser.add_argument("-r", "--remove-waters", action='store_true', help="Remove waters, salts etc.")

args = parser.parse_args()
print("analyse: ", args)

traj_in = args.trajectory
topol_in = args.protein
out_base = args.output

print('Reading trajectory', traj_in)
t = md.load(traj_in, top=topol_in)
t.image_molecules(inplace=True)

if args.remove_waters:
    print('Removing waters')
    t = t.atom_slice(t.top.select('not resname HOH POPC CL NA'))

print('Realigning')
prot = t.top.select('protein')
t.superpose(t[0], atom_indices=prot)

print('Writing re-imaged PDB', out_base + '.pdb')
t[0].save(out_base + '.pdb')

print('Writing re-imaged trajectory', out_base + '.dcd')
t.save(out_base + '.dcd')

topology = t.topology
# print(topology)
print('Number of frames:', t.n_frames)

# print('All residues: %s' % [residue for residue in t.topology.residues])

atoms = t.topology.select("chainid 1")
print(len(atoms), 'ligand atoms')
rmsds_lig = md.rmsd(t, t, frame=0, atom_indices=atoms, parallel=True, precentered=False)
# print(rmsds_lig)

atoms = t.topology.select("chainid 0 and backbone")
print(len(atoms), 'backbone atoms')
rmsds_bck = md.rmsd(t, t, frame=0, atom_indices=atoms, parallel=True, precentered=False)

fig = go.Figure()
fig.add_trace(go.Scatter(x=t.time, y=rmsds_lig, mode='lines', name='Ligand'))
fig.add_trace(go.Scatter(x=t.time, y=rmsds_bck, mode='lines', name='Backbone'))

fig.update_layout(title='Trajectory for ' + traj_in, xaxis_title='Frame', yaxis_title='RMSD')

file = out_base + '.svg'
print('Writing RMSD output to', file)
fig.write_image(file)