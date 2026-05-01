#!/usr/bin/env python3
"""Extract protein-only frames from a solvated OpenMM trajectory and compute CVs."""
import sys
import os
import numpy as np
import mdtraj as md

backup = sys.argv[1]  # e.g. steering_extended.bak.1150
out_dir = sys.argv[2]  # e.g. extracted_frames
os.makedirs(out_dir, exist_ok=True)

topo_pdb = f"{backup}/equilibrated.pdb"
nc_path = f"{backup}/production.nc"

print("Loading topology...")
topo = md.load(topo_pdb)
protein_sel = topo.topology.select("protein")
print(f"Protein atoms: {len(protein_sel)}, Total: {topo.n_atoms}")

print("Loading production trajectory...")
traj = md.load(nc_path, top=topo_pdb)
print(f"Loaded: {traj.n_frames} frames, {traj.n_atoms} atoms")

print("Extracting protein atoms...")
protein_traj = traj.atom_slice(protein_sel)
print(f"Protein: {protein_traj.n_frames} frames, {protein_traj.n_atoms} atoms")

# Compute CVs
print("Computing inter-domain CVs...")
top = protein_traj.topology
head_a = top.select("chainid 0 and resSeq 1 to 435")
head_b = top.select("chainid 1 and resSeq 1 to 352")
tail_a = top.select("chainid 0 and resSeq 780 to 9999")
tail_b = top.select("chainid 1 and resSeq 620 to 9999")

def domain_com(traj, indices):
    return traj.xyz[:, indices, :].mean(axis=1)

com_ha = domain_com(protein_traj, head_a)
com_hb = domain_com(protein_traj, head_b)
com_ta = domain_com(protein_traj, tail_a)
com_tb = domain_com(protein_traj, tail_b)

cv0 = np.linalg.norm(com_ha - com_tb, axis=1) * 10  # nm->A
cv1 = np.linalg.norm(com_hb - com_ta, axis=1) * 10
cv2 = np.linalg.norm(com_ha - com_hb, axis=1) * 10

cvs = np.column_stack([cv0, cv1, cv2])
np.save(f"{out_dir}/cvs.npy", cvs)

print(f"\nCV trajectory ({protein_traj.n_frames} frames):")
print(f"  CV0 (aH-bT): {cv0[0]:.1f} -> {cv0[-1]:.1f} A (max={cv0.max():.1f})")
print(f"  CV1 (bH-aT): {cv1[0]:.1f} -> {cv1[-1]:.1f} A (max={cv1.max():.1f})")
print(f"  CV2 (aH-bH): {cv2[0]:.1f} -> {cv2[-1]:.1f} A (max={cv2.max():.1f})")

# Save protein PDB frames
print(f"\nSaving {protein_traj.n_frames} protein PDB frames...")
for i in range(protein_traj.n_frames):
    protein_traj[i].save_pdb(f"{out_dir}/frame_{i:04d}.pdb")
    if (i + 1) % 100 == 0:
        print(f"  {i+1}/{protein_traj.n_frames}")

n_saved = len([f for f in os.listdir(out_dir) if f.endswith(".pdb")])
print(f"\nDone: {n_saved} frames, CVs in {out_dir}/cvs.npy")
sep = "=" * 60
print(f"\n{sep}")
print(f"CV0 range: {cv0.min():.1f} - {cv0.max():.1f} A")
print(f"CV1 range: {cv1.min():.1f} - {cv1.max():.1f} A")
print(f"Extension: CV0 max={cv0.max():.1f}A (target ~110A)")
last10 = cv0[-10:]
print(f"Still climbing? Last 10: {last10.mean():.1f} +/- {last10.std():.1f}")
print(sep)
