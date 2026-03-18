#!/usr/bin/env python3
"""Validate merged A5B1 tagged structure geometry.

Primary checks:
1) Expected chain count / required chains are present.
2) Tail proximity checks (A-tail->D and B-tail->C by default).
3) Optional chemistry checks at explicit attachment residues/atoms.
4) Optional exclusion checks to reject partner proximity to disallowed receptor regions.

This script reads a PDB file (no third-party deps) and exits non-zero on failure.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple


Coord = Tuple[float, float, float]
AtomRecord = Dict[str, object]


def parse_pdb_atoms(path: Path) -> Tuple[Dict[str, Dict[int, Coord]], List[AtomRecord]]:
    chains: Dict[str, Dict[int, Coord]] = {}
    atoms: List[AtomRecord] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        atom = line[12:16].strip()
        chain = line[21].strip()
        if not chain:
            continue
        try:
            residue = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except ValueError:
            continue
        xyz = (x, y, z)
        atoms.append(
            {
                "chain": chain,
                "residue": residue,
                "atom": atom,
                "xyz": xyz,
            }
        )
        if atom == "CA":
            chains.setdefault(chain, {})[residue] = xyz
    return chains, atoms


def euclidean(a: Coord, b: Coord) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def min_distance_to_chain(
    reference: Coord, target_residues: Dict[int, Coord]
) -> Tuple[Optional[float], Optional[int]]:
    best_distance: Optional[float] = None
    best_residue: Optional[int] = None
    for residue, xyz in target_residues.items():
        d = euclidean(reference, xyz)
        if best_distance is None or d < best_distance:
            best_distance = d
            best_residue = residue
    return best_distance, best_residue


def min_distance_between_atom_sets(
    source: Iterable[AtomRecord], target: Iterable[AtomRecord]
) -> Tuple[Optional[float], Optional[AtomRecord], Optional[AtomRecord]]:
    best_distance: Optional[float] = None
    best_source: Optional[AtomRecord] = None
    best_target: Optional[AtomRecord] = None
    target_list = list(target)
    for a in source:
        a_xyz = a.get("xyz")
        if not isinstance(a_xyz, tuple):
            continue
        for b in target_list:
            b_xyz = b.get("xyz")
            if not isinstance(b_xyz, tuple):
                continue
            d = euclidean(a_xyz, b_xyz)
            if best_distance is None or d < best_distance:
                best_distance = d
                best_source = a
                best_target = b
    return best_distance, best_source, best_target


def parse_residue_spec(spec: str) -> Set[int]:
    residues: Set[int] = set()
    raw = (spec or "").strip()
    if not raw:
        return residues
    for token in raw.split(","):
        item = token.strip()
        if not item:
            continue
        if "-" in item:
            left, right = item.split("-", 1)
            start = int(left.strip())
            end = int(right.strip())
            if end < start:
                start, end = end, start
            residues.update(range(start, end + 1))
        else:
            residues.add(int(item))
    return residues


def filter_atoms(
    atoms: List[AtomRecord],
    *,
    chain: Optional[str] = None,
    residue: Optional[int] = None,
    residues: Optional[Set[int]] = None,
    atom_names: Optional[Set[str]] = None,
) -> List[AtomRecord]:
    out: List[AtomRecord] = []
    for rec in atoms:
        if chain is not None and rec.get("chain") != chain:
            continue
        rec_res = rec.get("residue")
        if residue is not None and rec_res != residue:
            continue
        if residues is not None and rec_res not in residues:
            continue
        rec_atom = rec.get("atom")
        if atom_names is not None and rec_atom not in atom_names:
            continue
        out.append(rec)
    return out


def atom_ref(rec: Optional[AtomRecord]) -> Optional[str]:
    if rec is None:
        return None
    chain = rec.get("chain")
    residue = rec.get("residue")
    atom = rec.get("atom")
    return f"{chain}:{residue}:{atom}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", type=Path, required=True, help="Merged output PDB path.")
    parser.add_argument("--alpha-chain", default="A")
    parser.add_argument("--beta-chain", default="B")
    parser.add_argument(
        "--alpha-partner-chain",
        default="D",
        help="Partner chain expected to be attached near alpha tail (AviTag branch).",
    )
    parser.add_argument(
        "--beta-partner-chain",
        default="C",
        help="Partner chain expected to be attached near beta tail (SpyCatcher branch).",
    )
    parser.add_argument("--expected-chain-count", type=int, default=4)
    parser.add_argument(
        "--max-alpha-tail-distance",
        type=float,
        default=25.0,
        help="Maximum allowed CA distance (A tail -> alpha partner chain).",
    )
    parser.add_argument(
        "--max-beta-tail-distance",
        type=float,
        default=25.0,
        help="Maximum allowed CA distance (B tail -> beta partner chain).",
    )
    parser.add_argument(
        "--report-json",
        type=Path,
        default=None,
        help="Optional JSON report output path.",
    )
    parser.add_argument(
        "--alpha-attachment-residue",
        type=int,
        default=966,
        help="Alpha-chain residue expected to anchor AviTag/biotin branch.",
    )
    parser.add_argument(
        "--alpha-attachment-atom",
        default="NZ",
        help="Atom name on alpha attachment residue (default Lys NZ).",
    )
    parser.add_argument(
        "--max-alpha-attachment-distance",
        type=float,
        default=None,
        help="Optional max distance from alpha attachment atom to alpha partner chain (any atom).",
    )
    parser.add_argument(
        "--beta-attachment-residue",
        type=int,
        default=735,
        help="Beta-chain residue expected to react with SpyTag.",
    )
    parser.add_argument(
        "--beta-attachment-atom",
        default="NZ",
        help="Atom name on beta attachment residue (default Lys NZ).",
    )
    parser.add_argument(
        "--beta-partner-residue",
        type=int,
        default=10,
        help="Partner-chain residue on SpyTag expected to react with beta attachment residue.",
    )
    parser.add_argument(
        "--beta-partner-atoms",
        default="OD1,OD2,CG",
        help="Comma-separated atom names on beta partner residue for distance check.",
    )
    parser.add_argument(
        "--max-beta-attachment-distance",
        type=float,
        default=None,
        help="Optional max distance beta attachment atom -> beta partner reactive atoms.",
    )
    parser.add_argument(
        "--alpha-disallowed-residues",
        default="",
        help="Optional residue list/ranges on alpha chain (e.g. 780-820) that alpha partner should avoid.",
    )
    parser.add_argument(
        "--beta-disallowed-residues",
        default="",
        help="Optional residue list/ranges on beta chain that beta partner should avoid.",
    )
    parser.add_argument(
        "--min-alpha-partner-distance-to-alpha-disallowed",
        type=float,
        default=None,
        help="If set with --alpha-disallowed-residues, alpha partner must remain at least this far away.",
    )
    parser.add_argument(
        "--min-beta-partner-distance-to-beta-disallowed",
        type=float,
        default=None,
        help="If set with --beta-disallowed-residues, beta partner must remain at least this far away.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if not args.pdb.exists():
        raise FileNotFoundError(f"PDB not found: {args.pdb}")

    ca, atoms = parse_pdb_atoms(args.pdb)
    chain_ids = sorted(ca.keys())
    errors: List[str] = []

    if len(chain_ids) != args.expected_chain_count:
        errors.append(
            f"Expected {args.expected_chain_count} chains, found {len(chain_ids)} ({chain_ids})."
        )

    for required in (
        args.alpha_chain,
        args.beta_chain,
        args.alpha_partner_chain,
        args.beta_partner_chain,
    ):
        if required not in ca:
            errors.append(f"Required chain missing: {required}. Present: {chain_ids}")

    alpha_tail_res: Optional[int] = None
    beta_tail_res: Optional[int] = None
    alpha_tail_to_partner: Optional[float] = None
    beta_tail_to_partner: Optional[float] = None
    alpha_partner_res: Optional[int] = None
    beta_partner_res: Optional[int] = None
    alpha_attachment_distance: Optional[float] = None
    beta_attachment_distance: Optional[float] = None
    alpha_attachment_partner_atom: Optional[str] = None
    beta_attachment_partner_atom: Optional[str] = None
    alpha_disallowed_min_distance: Optional[float] = None
    beta_disallowed_min_distance: Optional[float] = None
    alpha_disallowed_partner_atom: Optional[str] = None
    alpha_disallowed_receptor_atom: Optional[str] = None
    beta_disallowed_partner_atom: Optional[str] = None
    beta_disallowed_receptor_atom: Optional[str] = None

    if args.alpha_chain in ca and args.alpha_partner_chain in ca:
        alpha_tail_res = max(ca[args.alpha_chain].keys())
        alpha_tail_xyz = ca[args.alpha_chain][alpha_tail_res]
        alpha_tail_to_partner, alpha_partner_res = min_distance_to_chain(
            alpha_tail_xyz, ca[args.alpha_partner_chain]
        )
        if (
            alpha_tail_to_partner is None
            or alpha_tail_to_partner > args.max_alpha_tail_distance
        ):
            errors.append(
                "Alpha tail distance check failed: "
                f"{alpha_tail_to_partner} A > {args.max_alpha_tail_distance} A"
            )

    if args.beta_chain in ca and args.beta_partner_chain in ca:
        beta_tail_res = max(ca[args.beta_chain].keys())
        beta_tail_xyz = ca[args.beta_chain][beta_tail_res]
        beta_tail_to_partner, beta_partner_res = min_distance_to_chain(
            beta_tail_xyz, ca[args.beta_partner_chain]
        )
        if beta_tail_to_partner is None or beta_tail_to_partner > args.max_beta_tail_distance:
            errors.append(
                "Beta tail distance check failed: "
                f"{beta_tail_to_partner} A > {args.max_beta_tail_distance} A"
            )

    # Optional explicit chemistry checks.
    alpha_attach_atom_set = filter_atoms(
        atoms,
        chain=args.alpha_chain,
        residue=args.alpha_attachment_residue,
        atom_names={args.alpha_attachment_atom},
    )
    alpha_partner_atoms = filter_atoms(atoms, chain=args.alpha_partner_chain)
    if args.max_alpha_attachment_distance is not None:
        if not alpha_attach_atom_set:
            errors.append(
                f"Missing alpha attachment atom {args.alpha_chain}:{args.alpha_attachment_residue}:{args.alpha_attachment_atom}"
            )
        elif not alpha_partner_atoms:
            errors.append(f"Missing alpha partner chain atoms: {args.alpha_partner_chain}")
        else:
            alpha_attachment_distance, _src, dst = min_distance_between_atom_sets(
                alpha_attach_atom_set, alpha_partner_atoms
            )
            alpha_attachment_partner_atom = atom_ref(dst)
            if (
                alpha_attachment_distance is None
                or alpha_attachment_distance > args.max_alpha_attachment_distance
            ):
                errors.append(
                    "Alpha attachment distance check failed: "
                    f"{alpha_attachment_distance} A > {args.max_alpha_attachment_distance} A"
                )

    beta_attach_atom_set = filter_atoms(
        atoms,
        chain=args.beta_chain,
        residue=args.beta_attachment_residue,
        atom_names={args.beta_attachment_atom},
    )
    beta_partner_atom_names = {
        x.strip() for x in args.beta_partner_atoms.split(",") if x.strip()
    }
    beta_partner_reactive_atoms = filter_atoms(
        atoms,
        chain=args.beta_partner_chain,
        residue=args.beta_partner_residue,
        atom_names=beta_partner_atom_names if beta_partner_atom_names else None,
    )
    if args.max_beta_attachment_distance is not None:
        if not beta_attach_atom_set:
            errors.append(
                f"Missing beta attachment atom {args.beta_chain}:{args.beta_attachment_residue}:{args.beta_attachment_atom}"
            )
        elif not beta_partner_reactive_atoms:
            errors.append(
                "Missing beta partner reactive atoms "
                f"{args.beta_partner_chain}:{args.beta_partner_residue}:{sorted(beta_partner_atom_names)}"
            )
        else:
            beta_attachment_distance, _src, dst = min_distance_between_atom_sets(
                beta_attach_atom_set, beta_partner_reactive_atoms
            )
            beta_attachment_partner_atom = atom_ref(dst)
            if (
                beta_attachment_distance is None
                or beta_attachment_distance > args.max_beta_attachment_distance
            ):
                errors.append(
                    "Beta attachment distance check failed: "
                    f"{beta_attachment_distance} A > {args.max_beta_attachment_distance} A"
                )

    # Optional exclusion checks against receptor regions.
    alpha_disallowed = parse_residue_spec(args.alpha_disallowed_residues)
    if args.min_alpha_partner_distance_to_alpha_disallowed is not None and alpha_disallowed:
        alpha_disallowed_atoms = filter_atoms(
            atoms,
            chain=args.alpha_chain,
            residues=alpha_disallowed,
        )
        if not alpha_partner_atoms:
            alpha_partner_atoms = filter_atoms(atoms, chain=args.alpha_partner_chain)
        if not alpha_disallowed_atoms:
            errors.append(
                f"No alpha disallowed residues found in chain {args.alpha_chain} for spec '{args.alpha_disallowed_residues}'"
            )
        elif not alpha_partner_atoms:
            errors.append(f"Missing alpha partner chain atoms: {args.alpha_partner_chain}")
        else:
            alpha_disallowed_min_distance, src, dst = min_distance_between_atom_sets(
                alpha_partner_atoms, alpha_disallowed_atoms
            )
            alpha_disallowed_partner_atom = atom_ref(src)
            alpha_disallowed_receptor_atom = atom_ref(dst)
            if (
                alpha_disallowed_min_distance is None
                or alpha_disallowed_min_distance
                < args.min_alpha_partner_distance_to_alpha_disallowed
            ):
                errors.append(
                    "Alpha disallowed-region check failed: "
                    f"{alpha_disallowed_min_distance} A < {args.min_alpha_partner_distance_to_alpha_disallowed} A"
                )

    beta_disallowed = parse_residue_spec(args.beta_disallowed_residues)
    if args.min_beta_partner_distance_to_beta_disallowed is not None and beta_disallowed:
        beta_disallowed_atoms = filter_atoms(
            atoms,
            chain=args.beta_chain,
            residues=beta_disallowed,
        )
        beta_partner_atoms = filter_atoms(atoms, chain=args.beta_partner_chain)
        if not beta_disallowed_atoms:
            errors.append(
                f"No beta disallowed residues found in chain {args.beta_chain} for spec '{args.beta_disallowed_residues}'"
            )
        elif not beta_partner_atoms:
            errors.append(f"Missing beta partner chain atoms: {args.beta_partner_chain}")
        else:
            beta_disallowed_min_distance, src, dst = min_distance_between_atom_sets(
                beta_partner_atoms, beta_disallowed_atoms
            )
            beta_disallowed_partner_atom = atom_ref(src)
            beta_disallowed_receptor_atom = atom_ref(dst)
            if (
                beta_disallowed_min_distance is None
                or beta_disallowed_min_distance
                < args.min_beta_partner_distance_to_beta_disallowed
            ):
                errors.append(
                    "Beta disallowed-region check failed: "
                    f"{beta_disallowed_min_distance} A < {args.min_beta_partner_distance_to_beta_disallowed} A"
                )

    report = {
        "pdb": str(args.pdb.resolve()),
        "chains_present": chain_ids,
        "expected_chain_count": args.expected_chain_count,
        "alpha_chain": args.alpha_chain,
        "beta_chain": args.beta_chain,
        "alpha_partner_chain": args.alpha_partner_chain,
        "beta_partner_chain": args.beta_partner_chain,
        "alpha_tail_residue": alpha_tail_res,
        "beta_tail_residue": beta_tail_res,
        "alpha_tail_to_partner_min_ca": alpha_tail_to_partner,
        "alpha_partner_nearest_residue": alpha_partner_res,
        "beta_tail_to_partner_min_ca": beta_tail_to_partner,
        "beta_partner_nearest_residue": beta_partner_res,
        "max_alpha_tail_distance": args.max_alpha_tail_distance,
        "max_beta_tail_distance": args.max_beta_tail_distance,
        "alpha_attachment_residue": args.alpha_attachment_residue,
        "alpha_attachment_atom": args.alpha_attachment_atom,
        "max_alpha_attachment_distance": args.max_alpha_attachment_distance,
        "alpha_attachment_distance": alpha_attachment_distance,
        "alpha_attachment_partner_atom": alpha_attachment_partner_atom,
        "beta_attachment_residue": args.beta_attachment_residue,
        "beta_attachment_atom": args.beta_attachment_atom,
        "beta_partner_residue": args.beta_partner_residue,
        "beta_partner_atoms": sorted(beta_partner_atom_names),
        "max_beta_attachment_distance": args.max_beta_attachment_distance,
        "beta_attachment_distance": beta_attachment_distance,
        "beta_attachment_partner_atom": beta_attachment_partner_atom,
        "alpha_disallowed_residues": sorted(alpha_disallowed),
        "min_alpha_partner_distance_to_alpha_disallowed": args.min_alpha_partner_distance_to_alpha_disallowed,
        "alpha_disallowed_min_distance": alpha_disallowed_min_distance,
        "alpha_disallowed_partner_atom": alpha_disallowed_partner_atom,
        "alpha_disallowed_receptor_atom": alpha_disallowed_receptor_atom,
        "beta_disallowed_residues": sorted(beta_disallowed),
        "min_beta_partner_distance_to_beta_disallowed": args.min_beta_partner_distance_to_beta_disallowed,
        "beta_disallowed_min_distance": beta_disallowed_min_distance,
        "beta_disallowed_partner_atom": beta_disallowed_partner_atom,
        "beta_disallowed_receptor_atom": beta_disallowed_receptor_atom,
        "pass": len(errors) == 0,
        "errors": errors,
    }

    if args.report_json is not None:
        args.report_json.parent.mkdir(parents=True, exist_ok=True)
        args.report_json.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print(json.dumps(report, indent=2))
    return 0 if report["pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
