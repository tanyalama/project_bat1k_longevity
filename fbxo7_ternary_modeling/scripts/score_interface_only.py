"""
Interface-only scoring of an ALREADY-RELAXED decoy PDB (preview while FastRelax finishes).
Mirrors the InterfaceAnalyzer half of relax_score_boltz_decoy.py exactly.
Usage: python score_interface_only.py <relaxed_pdb> <sample_name>
"""
import sys, csv
from pyrosetta import *
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.core.pose import append_pose_to_pose

init("-ex1 -ex2aro")

INPUT_PDB = sys.argv[1]
SAMPLE    = sys.argv[2]
WT_PDB    = "/work/pi_tlama_smith_edu/rosetta/results/WT_relaxed_best.pdb"
scorefxn  = ScoreFunctionFactory.create_score_function("ref2015")

def score_interface(pose, ca, cb):
    chains = pose.split_by_chain()
    pair = chains[ca].clone()
    append_pose_to_pose(pair, chains[cb], True)
    mover = InterfaceAnalyzerMover(1, False, scorefxn, True, False, False, False, True)
    mover.apply(pair)
    return {
        "dG_separated": mover.get_interface_dG(),
        "dSASA_int_A2": mover.get_interface_delta_sasa(),
        "delta_unsatHbonds": mover.get_interface_delta_hbond_unsat(),
        "packstat": mover.get_interface_packstat(),
        "nres_int": mover.get_num_interface_residues(),
    }

interfaces = {"A_B": ("FBXO7-PINK1", 1, 2), "A_C": ("FBXO7-PSMF1", 1, 3)}
wt_pose  = pose_from_file(WT_PDB)
mut_pose = pose_from_file(INPUT_PDB)
print(f"loaded {INPUT_PDB}: res109 = {mut_pose.residue(109).name3()}", flush=True)

rows = []
for key, (label, ca, cb) in interfaces.items():
    wt = score_interface(wt_pose, ca, cb)
    rows.append({"sample": "WT", "interface": key, "label": label,
                 "dG_separated": wt["dG_separated"], "ddG": 0.0,
                 "dSASA_int_A2": wt["dSASA_int_A2"], "delta_unsatHbonds": wt["delta_unsatHbonds"],
                 "packstat": wt["packstat"], "nres_int": wt["nres_int"], "pdb_file": WT_PDB})
for key, (label, ca, cb) in interfaces.items():
    wt = score_interface(wt_pose, ca, cb)
    mut = score_interface(mut_pose, ca, cb)
    ddg = mut["dG_separated"] - wt["dG_separated"]
    rows.append({"sample": SAMPLE, "interface": key, "label": label,
                 "dG_separated": mut["dG_separated"], "ddG": ddg,
                 "dSASA_int_A2": mut["dSASA_int_A2"], "delta_unsatHbonds": mut["delta_unsatHbonds"],
                 "packstat": mut["packstat"], "nres_int": mut["nres_int"], "pdb_file": INPUT_PDB})
    print(f"{SAMPLE}\t{key}\t{label}\tdG={mut['dG_separated']:.3f}\tddG={ddg:.3f}", flush=True)

with open(f"./{SAMPLE}_preview_interface_metrics.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["sample","interface","label","dG_separated","ddG",
        "dSASA_int_A2","delta_unsatHbonds","packstat","nres_int","pdb_file"])
    w.writeheader(); w.writerows(rows)
print("DONE", SAMPLE, flush=True)
