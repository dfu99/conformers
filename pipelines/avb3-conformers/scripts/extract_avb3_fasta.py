#!/usr/bin/env python3
"""Extract AVB3 chain sequences from a PDB and write FASTA for AF2."""
import sys
from pathlib import Path

def main():
    if len(sys.argv) < 3:
        print("Usage: extract_avb3_fasta.py <pdb_path> <fasta_dir>")
        sys.exit(1)

    pdb_path = Path(sys.argv[1])
    fasta_dir = Path(sys.argv[2])
    fasta_dir.mkdir(parents=True, exist_ok=True)

    aa3to1 = {
        "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
        "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
        "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M",
        "HSD":"H","HSE":"H","HSP":"H",
    }
    chains = {}
    seen = set()
    for line in pdb_path.read_text().splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        cid = line[21]
        resseq = int(line[22:26])
        key = (cid, resseq)
        if key in seen:
            continue
        seen.add(key)
        chains.setdefault(cid, []).append((resseq, aa3to1.get(line[17:20].strip(), "X")))

    seqs = {c: "".join(aa for _, aa in sorted(res)) for c, res in sorted(chains.items())}
    fasta_path = fasta_dir / "avb3_multimer.fasta"
    with open(fasta_path, "w") as f:
        for c, s in sorted(seqs.items()):
            f.write(f">chain_{c}\n{s}\n")
    print(f"Wrote {fasta_path}")
    for c, s in sorted(seqs.items()):
        print(f"  Chain {c}: {len(s)} residues")

if __name__ == "__main__":
    main()
