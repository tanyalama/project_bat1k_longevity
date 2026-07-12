"""
Constrained FastRelax + InterfaceAnalyzer scoring for a Boltz-predicted mutant decoy.
Mirrors the established smith pipeline (05_test_s109p_5seed.py + 06_batch_scoring.py)
but SKIPS MutateResidue because the Boltz decoy already carries the mutation.

Usage: python relax_score_boltz_decoy.py <input_pdb> <sample_name>
"""
import sys, os, csv
from pyrosetta import *
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.constraint_generator import (
    CoordinateConstraintGenerator, AddConstraints)
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.core.pose import append_pose_to_pose
import pyrosetta.rosetta.core.scoring as scoring

init("-ex1 -ex2aro")

INPUT_PDB = sys.argv[1]
SAMPLE    = sys.argv[2]
WT_PDB    = "/work/pi_tlama_smith_edu/rosetta/results/WT_relaxed_best.pdb"
SD        = 0.5
N_SEEDS   = 5

scorefxn = ScoreFunctionFactory.create_score_function("ref2015")

# ---------- constrained FastRelax, keep lowest-energy decoy ----------
best_score = None
best_pose  = None
relax_rows = []
for seed in range(1, N_SEEDS + 1):
    pose = pose_from_file(INPUT_PDB)
    print(f"\n[seed {seed}] res109 = {pose.residue(109).name3()}  loaded {pose.total_residue()} res", flush=True)

    cst_gen = CoordinateConstraintGenerator()
    cst_gen.set_sd(SD)
    cst_gen.set_ca_only(True)
    add = AddConstraints()
    add.add_generator(cst_gen)
    add.apply(pose)

    sf = ScoreFunctionFactory.create_score_function("ref2015")
    sf.set_weight(scoring.coordinate_constraint, 1.0)
    init_s = sf(pose)

    relax = FastRelax()
    relax.set_scorefxn(sf)
    relax.apply(pose)
    final_s = sf(pose)
    print(f"[seed {seed}] init={init_s:.2f} final={final_s:.2f}", flush=True)

    out_pdb = f"./{SAMPLE}_relaxed_CA_sd{SD}_seed{seed}.pdb"
    pose.dump_pdb(out_pdb)
    relax_rows.append({
        "sample": SAMPLE, "seed": seed, "constraint_sd": SD,
        "initial_score_with_constraints": init_s,
        "final_score_with_constraints": final_s,
        "score_change": final_s - init_s, "output_pdb": out_pdb,
    })
    if best_score is None or final_s < best_score:
        best_score = final_s
        best_pose  = pose.clone()
        best_seed  = seed

best_pdb = f"./{SAMPLE}_relaxed_best.pdb"
best_pose.dump_pdb(best_pdb)
print(f"\nBest decoy: seed {best_seed}, constrained score {best_score:.2f} -> {best_pdb}", flush=True)

with open(f"./{SAMPLE}_relax_summary.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(relax_rows[0].keys()))
    w.writeheader(); w.writerows(relax_rows)

# ---------- interface scoring vs WT ----------
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
wt_pose = pose_from_file(WT_PDB)

rows = []
for key, (label, ca, cb) in interfaces.items():
    wt = score_interface(wt_pose, ca, cb)
    rows.append({"sample": "WT", "interface": key, "label": label,
                 "dG_separated": wt["dG_separated"], "ddG": 0.0,
                 "dSASA_int_A2": wt["dSASA_int_A2"], "delta_unsatHbonds": wt["delta_unsatHbonds"],
                 "packstat": wt["packstat"], "nres_int": wt["nres_int"], "pdb_file": WT_PDB})
for key, (label, ca, cb) in interfaces.items():
    wt = score_interface(wt_pose, ca, cb)
    mut = score_interface(best_pose, ca, cb)
    ddg = mut["dG_separated"] - wt["dG_separated"]
    rows.append({"sample": SAMPLE, "interface": key, "label": label,
                 "dG_separated": mut["dG_separated"], "ddG": ddg,
                 "dSASA_int_A2": mut["dSASA_int_A2"], "delta_unsatHbonds": mut["delta_unsatHbonds"],
                 "packstat": mut["packstat"], "nres_int": mut["nres_int"], "pdb_file": best_pdb})
    print(f"{SAMPLE}\t{key}\t{label}\tdG={mut['dG_separated']:.3f}\tddG={ddg:.3f}", flush=True)

with open(f"./{SAMPLE}_interface_metrics.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["sample","interface","label","dG_separated","ddG",
        "dSASA_int_A2","delta_unsatHbonds","packstat","nres_int","pdb_file"])
    w.writeheader(); w.writerows(rows)
print("DONE", SAMPLE, flush=True)
