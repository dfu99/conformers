#!/usr/bin/env python3
"""Register a custom Protenix template entry without running inference.

This performs only registration steps:
1) Copy template mmCIF into Protenix mmcif DB layout: mmcif/<id[1:3]>/<id>.cif.gz
2) Add/overwrite template entry release date in common/release_date_cache.json
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import shutil
from pathlib import Path

ENTRY_ID_RE = re.compile(r"^[A-Za-z0-9_]{3,}$")
DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")


def resolve_mmcif_dir(explicit: Path | None) -> Path:
    if explicit is not None:
        return explicit.resolve()
    root = os.environ.get("PROTENIX_ROOT_DIR", "").strip()
    if not root:
        raise RuntimeError(
            "Cannot resolve mmcif dir. Set --template_mmcif_dir or PROTENIX_ROOT_DIR."
        )
    return (Path(root) / "mmcif").resolve()


def resolve_release_dates_path(explicit: Path | None) -> Path:
    if explicit is not None:
        return explicit.resolve()
    root = os.environ.get("PROTENIX_ROOT_DIR", "").strip()
    if not root:
        raise RuntimeError(
            "Cannot resolve release_date_cache.json. "
            "Set --release_dates_path or PROTENIX_ROOT_DIR."
        )
    return (Path(root) / "common" / "release_date_cache.json").resolve()


def register_mmcif(template_cif: Path, entry_id: str, mmcif_dir: Path) -> Path:
    if not template_cif.exists():
        raise FileNotFoundError(f"Template CIF not found: {template_cif}")

    subdir = mmcif_dir / entry_id[1:3]
    subdir.mkdir(parents=True, exist_ok=True)
    target = subdir / f"{entry_id}.cif.gz"

    with template_cif.open("rb") as src, gzip.open(target, "wb") as dst:
        shutil.copyfileobj(src, dst)

    return target.resolve()


def register_release_date(entry_id: str, release_date: str, release_dates_path: Path) -> Path:
    if not DATE_RE.fullmatch(release_date):
        raise ValueError("--release_date must have format YYYY-MM-DD")

    release_dates_path.parent.mkdir(parents=True, exist_ok=True)
    if release_dates_path.exists():
        data = json.loads(release_dates_path.read_text(encoding="utf-8"))
        if not isinstance(data, dict):
            raise ValueError(f"Expected JSON object in {release_dates_path}")
    else:
        data = {}

    data[entry_id.lower()] = release_date
    release_dates_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    return release_dates_path.resolve()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Register custom Protenix template files (mmcif + release date cache)."
    )
    parser.add_argument(
        "--template_cif",
        type=Path,
        required=True,
        help="Path to template .cif file to register.",
    )
    parser.add_argument(
        "--template_entry_id",
        type=str,
        default="s090",
        help="Template entry ID (e.g. s090).",
    )
    parser.add_argument(
        "--template_mmcif_dir",
        type=Path,
        default=None,
        help="Protenix mmcif DB root (default: $PROTENIX_ROOT_DIR/mmcif).",
    )
    parser.add_argument(
        "--release_dates_path",
        type=Path,
        default=None,
        help=(
            "Path to release_date_cache.json "
            "(default: $PROTENIX_ROOT_DIR/common/release_date_cache.json)."
        ),
    )
    parser.add_argument(
        "--release_date",
        type=str,
        default="1999-01-01",
        help="Release date to assign for this custom template entry.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    entry_id = args.template_entry_id.strip().lower()
    if not ENTRY_ID_RE.fullmatch(entry_id):
        raise ValueError("--template_entry_id must match [A-Za-z0-9_]{3,}")

    template_cif = args.template_cif.resolve()
    mmcif_dir = resolve_mmcif_dir(args.template_mmcif_dir)
    release_dates_path = resolve_release_dates_path(args.release_dates_path)

    mmcif_target = register_mmcif(template_cif=template_cif, entry_id=entry_id, mmcif_dir=mmcif_dir)
    release_dates_target = register_release_date(
        entry_id=entry_id,
        release_date=args.release_date,
        release_dates_path=release_dates_path,
    )

    print("Registered template entry:")
    print(f"  entry_id:            {entry_id}")
    print(f"  template_cif:        {template_cif}")
    print(f"  mmcif_registered_as: {mmcif_target}")
    print(f"  release_dates_path:  {release_dates_target}")
    print(f"  release_date:        {args.release_date}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
