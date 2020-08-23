import sys
import mdtraj as md
import plotly.graph_objects as go
import plotly.express as px

if len(sys.argv) != 4:
    print('Usage: python analyse.py trajectory.dcd topology.pdb output')
    print('Analyse MD trajectory')
    exit(1)

traj_in = sys.argv[1]
topol_in = sys.argv[2]
outfile = sys.argv[3]

print('Reading trajectory', traj_in)
t = md.load(traj_in, top=topol_in)

topology = t.topology
print(topology)
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

file = outfile + '.svg'
print('Writing output to', file)
fig.write_image(file)