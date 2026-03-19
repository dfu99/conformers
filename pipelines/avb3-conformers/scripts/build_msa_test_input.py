#!/usr/bin/env python3
"""Build a single Protenix input JSON with pre-computed MSA paths."""
from __future__ import annotations
import argparse, json
from pathlib import Path

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ref-pdb", type=Path, required=True)
    p.add_argument("--msa-a", type=str, required=True)
    p.add_argument("--msa-b", type=str, required=True)
    p.add_argument("--label", type=str, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()

    aa3to1 = {
        "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
        "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
        "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M",
        "HSD":"H","HSE":"H","HSP":"H",
    }
    chains: dict[str, list] = {}
    seen = set()
    for line in args.ref_pdb.read_text().splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        cid = line[21]
        resseq = int(line[22:26])
        key = (cid, resseq)
        if key in seen:
            continue
        seen.add(key)
        chains.setdefault(cid, []).append((resseq, aa3to1.get(line[17:20].strip(), "X")))

    chain_ids = sorted(chains.keys())
    seqs = {c: "".join(aa for _, aa in sorted(res)) for c, res in chains.items()}

    msa_map = {chain_ids[0]: args.msa_a, chain_ids[1]: args.msa_b}

    payload = [{"name": f"avb3_msa_{args.label}", "covalent_bonds": [], "sequences": [
        {"proteinChain": {"count": 1, "sequence": seqs[c], "modifications": [],
         "msa": {"precomputed_msa_dir": "", "pairedMsaPath": "", "unpairedMsaPath": msa_map[c]}}}
        for c in chain_ids
    ]}]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2))
    for c in chain_ids:
        print(f"  Chain {c}: {len(seqs[c])} residues, MSA: {msa_map[c]}")

if __name__ == "__main__":
    main()
