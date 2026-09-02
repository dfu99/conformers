#!/usr/bin/env python3
"""Discover + VERIFY full-ectodomain integrin heterodimer structures in the PDB.

Motivation (tasks/lessons.md, obj-070): a hand-quoted PDB ID ("6WOV = alpha5beta1
bent cryo-EM") turned out to be the ryanodine receptor. Never trust a recalled
accession. This script therefore does not hardcode a candidate list from memory --
it queries RCSB for entities that carry an integrin BETA-subunit UniProt accession,
then verifies every hit against the deposited coordinates themselves.

A structure qualifies as a morphable bent ectodomain only if, in the actual file:
  * exactly one alpha-like and one beta-like protein chain are identifiable,
  * the alpha chain has >= MIN_ALPHA_CA observed CA atoms (needs the legs:
    beta-propeller + thigh + calf-1 + calf-2, ~900 res in a full ectodomain),
  * the beta chain has >= MIN_BETA_CA observed CA atoms (betaI + hybrid + PSI +
    I-EGF1-4 + beta-tail, ~600 res),
  * the chain is long enough end-to-end to actually possess a genu to swing.

Headpiece-only structures (the majority of integrin entries) are reported but
excluded -- with no legs there is nothing to extend.

Pure stdlib + numpy. Network required (RCSB Search + Data API, files.rcsb.org).
"""
import argparse
import json
import os
import time
import urllib.request
import urllib.error

import numpy as np

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{}"
ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/{}/{}"
PDB_URL = "https://files.rcsb.org/download/{}.pdb"

# Human integrin subunit UniProt accessions.
BETA_ACC = {
    "P05556": "beta1", "P05107": "beta2", "P05106": "beta3", "P18084": "beta5",
    "P18564": "beta6", "P26010": "beta7", "P26012": "beta8", "P16144": "beta4",
}
ALPHA_ACC = {
    "P06756": "alphaV", "P08514": "alphaIIb", "P08648": "alpha5", "P20701": "alphaL",
    "P20702": "alphaX", "P11215": "alphaM", "P13612": "alpha4", "P56199": "alpha1",
    "P17301": "alpha2", "P26006": "alpha3", "P23229": "alpha6", "Q13683": "alpha7",
    "P53708": "alpha8", "Q13797": "alpha9", "O75578": "alpha10", "Q9UKX5": "alpha11",
    "P26006": "alpha3", "P13612": "alpha4", "P38570": "alphaE", "P08648": "alpha5",
}

# A full ectodomain must show most of both chains. Bent crystals have disordered
# loops, so these floors sit below the nominal lengths (alpha ~960, beta ~690).
MIN_ALPHA_CA = 780
MIN_BETA_CA = 520
MIN_LONG_AXIS_A = 90.0  # a headpiece alone spans ~70-80 A; legged forms exceed this


def post_json(url, payload, retries=3):
    data = json.dumps(payload).encode()
    req = urllib.request.Request(url, data=data,
                                 headers={"Content-Type": "application/json"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=60) as fh:
                return json.load(fh)
        except urllib.error.HTTPError as e:
            if e.code == 204:
                return {"result_set": []}
            if attempt == retries - 1:
                raise
        except Exception:
            if attempt == retries - 1:
                raise
        time.sleep(2 * (attempt + 1))


def get_json(url, retries=3):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                return json.load(fh)
        except Exception:
            if attempt == retries - 1:
                return None
            time.sleep(2 * (attempt + 1))
    return None


def search_beta_entities():
    """Return entry IDs having a polymer entity mapped to an integrin beta accession."""
    payload = {
        "query": {
            "type": "terminal", "service": "text",
            "parameters": {
                "attribute": ("rcsb_polymer_entity_container_identifiers"
                              ".reference_sequence_identifiers.database_accession"),
                "operator": "in",
                "value": sorted(BETA_ACC),
            },
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": 2000},
            "results_verbosity": "compact",
        },
    }
    res = post_json(SEARCH_URL, payload)
    return list(res.get("result_set", []))


def entry_entities(pdb_id):
    """Return [{entity_id, name, accessions, auth_chains, sample_len}] for an entry."""
    entry = get_json(ENTRY_URL.format(pdb_id))
    if not entry:
        return None, None
    title = entry.get("struct", {}).get("title", "")
    ids = entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])
    out = []
    for eid in ids:
        ent = get_json(ENTITY_URL.format(pdb_id, eid))
        if not ent:
            continue
        cont = ent.get("rcsb_polymer_entity_container_identifiers", {})
        refs = cont.get("reference_sequence_identifiers") or []
        accs = [r.get("database_accession") for r in refs if r.get("database_accession")]
        out.append({
            "entity_id": eid,
            "name": (ent.get("rcsb_polymer_entity", {}) or {}).get("pdbx_description", ""),
            "accessions": accs,
            "auth_chains": cont.get("auth_asym_ids", []),
            "sample_len": (ent.get("entity_poly", {}) or {}).get(
                "rcsb_sample_sequence_length"),
        })
    return title, out


def classify(entity):
    """Return ('alpha'|'beta'|None, subunit_label) for a polymer entity."""
    for acc in entity["accessions"]:
        if acc in ALPHA_ACC:
            return "alpha", ALPHA_ACC[acc]
        if acc in BETA_ACC:
            return "beta", BETA_ACC[acc]
    nm = (entity["name"] or "").lower()
    if "integrin" in nm and "alpha" in nm:
        return "alpha", entity["name"]
    if "integrin" in nm and "beta" in nm:
        return "beta", entity["name"]
    return None, None


def download_pdb(pdb_id, dest_dir):
    os.makedirs(dest_dir, exist_ok=True)
    path = os.path.join(dest_dir, f"{pdb_id}.pdb")
    if os.path.exists(path) and os.path.getsize(path) > 10000:
        return path
    try:
        urllib.request.urlretrieve(PDB_URL.format(pdb_id), path)
    except Exception:
        return None
    return path if os.path.exists(path) else None


def ca_by_chain(path):
    """{auth_chain: [(resseq, xyz)]} for first-altloc CA atoms of the first model."""
    chains = {}
    seen = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            if line[16] not in (" ", "A"):
                continue
            ch = line[21]
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            if (ch, resseq) in seen:
                continue
            seen.add((ch, resseq))
            chains.setdefault(ch, []).append(
                (resseq, np.array([float(line[30:38]), float(line[38:46]),
                                   float(line[46:54])])))
    for ch in chains:
        chains[ch].sort(key=lambda t: t[0])
    return chains


def long_axis_extent(pts):
    c = pts.mean(axis=0)
    _, _, vt = np.linalg.svd(pts - c, full_matrices=False)
    proj = (pts - c) @ vt[0]
    return float(proj.max() - proj.min())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb-dir", default="inputs/variants")
    ap.add_argument("--out-json", default="results/variants/discovered_ectodomains.json")
    ap.add_argument("--max-entries", type=int, default=0, help="0 = no cap")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_json), exist_ok=True)

    print("[discover] querying RCSB for integrin beta-subunit entities...", flush=True)
    entries = search_beta_entities()
    if args.max_entries:
        entries = entries[: args.max_entries]
    print(f"[discover] {len(entries)} candidate entries", flush=True)

    accepted, rejected = [], []
    for i, pdb_id in enumerate(entries, 1):
        title, ents = entry_entities(pdb_id)
        if not ents:
            rejected.append({"pdb_id": pdb_id, "reason": "metadata fetch failed"})
            continue
        alphas = [e for e in ents if classify(e)[0] == "alpha"]
        betas = [e for e in ents if classify(e)[0] == "beta"]
        if not alphas or not betas:
            rejected.append({"pdb_id": pdb_id, "title": title,
                             "reason": "no alpha+beta integrin pair"})
            continue
        # cheap pre-filter on construct length before spending a download
        a_len = max([e["sample_len"] or 0 for e in alphas])
        b_len = max([e["sample_len"] or 0 for e in betas])
        if a_len < MIN_ALPHA_CA or b_len < MIN_BETA_CA:
            rejected.append({"pdb_id": pdb_id, "title": title,
                             "reason": f"construct too short (alpha {a_len}, beta {b_len})"})
            continue

        path = download_pdb(pdb_id, args.pdb_dir)
        if not path:
            rejected.append({"pdb_id": pdb_id, "title": title,
                             "reason": "PDB-format download unavailable (likely large/cryoEM-only mmCIF)"})
            continue
        chains = ca_by_chain(path)

        def best_chain(ent_list):
            cands = []
            for e in ent_list:
                for ch in e["auth_chains"]:
                    if ch in chains:
                        cands.append((len(chains[ch]), ch, e))
            return max(cands) if cands else None

        a_best, b_best = best_chain(alphas), best_chain(betas)
        if not a_best or not b_best:
            rejected.append({"pdb_id": pdb_id, "title": title,
                             "reason": "entity chains absent from coordinates"})
            continue
        a_n, a_ch, a_ent = a_best
        b_n, b_ch, b_ent = b_best
        if a_n < MIN_ALPHA_CA or b_n < MIN_BETA_CA:
            rejected.append({"pdb_id": pdb_id, "title": title,
                             "reason": f"headpiece-only / partial (observed CA: alpha {a_n}, beta {b_n})"})
            continue

        pts = np.array([xyz for ch in (a_ch, b_ch) for _, xyz in chains[ch]])
        extent = long_axis_extent(pts)
        if extent < MIN_LONG_AXIS_A:
            rejected.append({"pdb_id": pdb_id, "title": title,
                             "reason": f"long axis {extent:.0f} A < {MIN_LONG_AXIS_A} (no legs)"})
            continue

        a_res = [r for r, _ in chains[a_ch]]
        b_res = [r for r, _ in chains[b_ch]]
        rec = {
            "pdb_id": pdb_id, "title": title,
            "alpha": {"chain": a_ch, "subunit": classify(a_ent)[1], "name": a_ent["name"],
                      "accessions": a_ent["accessions"], "observed_ca": a_n,
                      "resseq_min": a_res[0], "resseq_max": a_res[-1]},
            "beta": {"chain": b_ch, "subunit": classify(b_ent)[1], "name": b_ent["name"],
                     "accessions": b_ent["accessions"], "observed_ca": b_n,
                     "resseq_min": b_res[0], "resseq_max": b_res[-1]},
            "long_axis_extent_A": round(extent, 1),
            "pdb_path": path,
        }
        rec["variant"] = f"{rec['alpha']['subunit']}{rec['beta']['subunit']}"
        accepted.append(rec)
        print(f"[discover] ACCEPT {pdb_id} {rec['variant']:<16} "
              f"alphaCA={a_n} betaCA={b_n} extent={extent:.0f} A  | {title[:60]}",
              flush=True)
        if i % 25 == 0:
            print(f"[discover] ...{i}/{len(entries)} scanned", flush=True)

    # one representative per variant: prefer the most complete (most observed CA)
    by_variant = {}
    for r in accepted:
        v = r["variant"]
        score = r["alpha"]["observed_ca"] + r["beta"]["observed_ca"]
        if v not in by_variant or score > by_variant[v]["_score"]:
            by_variant[v] = {**r, "_score": score}

    result = {
        "n_candidates": len(entries),
        "n_accepted": len(accepted),
        "criteria": {"min_alpha_ca": MIN_ALPHA_CA, "min_beta_ca": MIN_BETA_CA,
                     "min_long_axis_A": MIN_LONG_AXIS_A},
        "accepted": accepted,
        "representatives": {k: {kk: vv for kk, vv in v.items() if kk != "_score"}
                            for k, v in by_variant.items()},
        "rejected_sample": rejected[:60],
        "n_rejected": len(rejected),
    }
    with open(args.out_json, "w") as fh:
        json.dump(result, fh, indent=2)

    print(f"\n[discover] {len(accepted)} full-ectodomain structures across "
          f"{len(by_variant)} variants")
    for v, r in sorted(by_variant.items()):
        print(f"    {v:<18} {r['pdb_id']}  alphaCA={r['alpha']['observed_ca']:<5}"
              f" betaCA={r['beta']['observed_ca']:<5} extent={r['long_axis_extent_A']} A")
    print(f"[discover] wrote {args.out_json}")


if __name__ == "__main__":
    main()
