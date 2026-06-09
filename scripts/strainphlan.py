# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
summarize_strainphlan.py  (v4)
Summarizes StrainPhlAn4 results across all clade output directories,
optionally joining species/genus names from output_sgb_names.tsv.

Usage:
    python summarize_strainphlan.py \\
        --root   /path/to/Output \\
        --names  /path/to/output_sgb_names.tsv \\
        --out    strainphlan_summary.tsv

Columns produced:
    clade
    species                     e.g. "Bacteroides uniformis"
    genus                       e.g. "Bacteroides"
    full_lineage                full d__;p__;...;s__ string
    n_samples_initial
    n_samples_in_tree
    n_samples_excluded
    n_markers_available
    n_markers_used
    n_references
    n_samples_polymorphic_data
    mean_pct_polymorphic
    max_pct_polymorphic
    sample_max_pct_polymorphic
"""

import re
import argparse
import csv
from pathlib import Path


# ---------------------------------------------------------------------------
# Load SGB name table
# ---------------------------------------------------------------------------

def load_sgb_names(names_path):
    """
    Returns dict: clade_id -> {"species": ..., "genus": ..., "full_lineage": ...}
    Handles full lineage strings like:
        d__Bacteria;p__Bacteroidota;...;g__Bacteroides;s__Bacteroides uniformis
    """
    names = {}
    if names_path is None or not Path(names_path).exists():
        return names

    with open(names_path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            clade = row.get("Clade", "").strip()
            lineage = row.get("Name", "").strip()
            if not clade:
                continue

            species = ""
            genus   = ""
            for part in lineage.split(";"):
                part = part.strip()
                if part.startswith("s__"):
                    species = part[3:]
                elif part.startswith("g__"):
                    genus = part[3:]

            names[clade] = {
                "species":      species,
                "genus":        genus,
                "full_lineage": lineage,
            }
    return names


# ---------------------------------------------------------------------------
# .info parsing
# ---------------------------------------------------------------------------

def parse_info_file(info_path):
    result = {
        "n_samples_initial":   None,
        "n_samples_in_tree":   None,
        "n_markers_available": None,
        "n_markers_used":      None,
        "n_references":        None,
    }
    if not info_path.exists():
        return result

    text = info_path.read_text(errors="replace")

    patterns = {
        "n_samples_initial":   r"Number of samples:\s*(\d+)",
        "n_samples_in_tree":   r"Number of samples after filtering:\s*(\d+)",
        "n_markers_available": r"Number of available markers for the clade:\s*(\d+)",
        "n_markers_used":      r"Number of markers selected after filtering:\s*(\d+)",
        "n_references":        r"Number of references after filtering:\s*(\d+)",
    }

    for key, pat in patterns.items():
        m = re.search(pat, text, re.I)
        if m:
            result[key] = int(m.group(1))

    return result


# ---------------------------------------------------------------------------
# .polymorphic parsing
# ---------------------------------------------------------------------------

def parse_polymorphic_file(poly_path):
    result = {
        "n_samples_polymorphic_data":  0,
        "mean_pct_polymorphic":        None,
        "max_pct_polymorphic":         None,
        "sample_max_pct_polymorphic":  None,
    }
    if not poly_path.exists():
        return result

    text = poly_path.read_text(errors="replace").strip()
    if not text:
        return result

    lines = text.splitlines()

    header_idx = None
    for i, line in enumerate(lines):
        fields = line.strip().split("\t")
        if fields[0].lower() == "sample" or "percentage_of_polymorphic_sites" in line:
            header_idx = i
            break

    if header_idx is None:
        result["n_samples_polymorphic_data"] = sum(1 for l in lines if l.strip())
        return result

    headers    = lines[header_idx].strip().split("\t")
    data_lines = [l for l in lines[header_idx + 1:] if l.strip()]

    try:
        pct_col    = headers.index("percentage_of_polymorphic_sites")
        sample_col = headers.index("sample")
    except ValueError:
        result["n_samples_polymorphic_data"] = len(data_lines)
        return result

    pct_values = []
    samples    = []
    for line in data_lines:
        fields = line.strip().split("\t")
        if len(fields) <= max(pct_col, sample_col):
            continue
        try:
            pct = float(fields[pct_col])
            pct_values.append(pct)
            samples.append(fields[sample_col])
        except ValueError:
            continue

    result["n_samples_polymorphic_data"] = len(pct_values)
    if pct_values:
        result["mean_pct_polymorphic"] = round(sum(pct_values) / len(pct_values), 4)
        max_pct = max(pct_values)
        result["max_pct_polymorphic"]        = round(max_pct, 4)
        result["sample_max_pct_polymorphic"] = samples[pct_values.index(max_pct)]

    return result


# ---------------------------------------------------------------------------
# Tree tip counting
# ---------------------------------------------------------------------------

def count_tree_tips(tree_path):
    if not tree_path.exists():
        return None
    text = tree_path.read_text(errors="replace").strip()
    leaves = re.findall(r"[,(]([^,();:]+)(?::[0-9eE.\-]+)?", text)
    all_leaves = {x.strip() for x in leaves
                  if x.strip() and not re.fullmatch(r"[\d.eE\-]+", x.strip())}
    return len(all_leaves) if all_leaves else None


# ---------------------------------------------------------------------------
# Per-clade summary
# ---------------------------------------------------------------------------

def summarize_clade(clade_dir, sgb_names):
    clade_name = clade_dir.name
    info_path  = clade_dir / f"{clade_name}.info"
    poly_path  = clade_dir / f"{clade_name}.polymorphic"
    tree_path  = clade_dir / f"RAxML_bestTree.{clade_name}.StrainPhlAn4.tre"

    info  = parse_info_file(info_path)
    poly  = parse_polymorphic_file(poly_path)
    name  = sgb_names.get(clade_name, {"species": "", "genus": "", "full_lineage": ""})

    n_tree_tips = count_tree_tips(tree_path)
    n_in_tree   = n_tree_tips if n_tree_tips is not None else info["n_samples_in_tree"]

    n_excluded = None
    if info["n_samples_initial"] is not None and n_in_tree is not None:
        n_excluded = info["n_samples_initial"] - n_in_tree

    return {
        "clade":                        clade_name,
        "species":                      name["species"],
        "genus":                        name["genus"],
        "full_lineage":                 name["full_lineage"],
        "n_samples_initial":            info["n_samples_initial"],
        "n_samples_in_tree":            n_in_tree,
        "n_samples_excluded":           n_excluded,
        "n_markers_available":          info["n_markers_available"],
        "n_markers_used":               info["n_markers_used"],
        "n_references":                 info["n_references"],
        "n_samples_polymorphic_data":   poly["n_samples_polymorphic_data"],
        "mean_pct_polymorphic":         poly["mean_pct_polymorphic"],
        "max_pct_polymorphic":          poly["max_pct_polymorphic"],
        "sample_max_pct_polymorphic":   poly["sample_max_pct_polymorphic"],
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

FIELDNAMES = [
    "clade",
    "species",
    "genus",
    "full_lineage",
    "n_samples_initial",
    "n_samples_in_tree",
    "n_samples_excluded",
    "n_markers_available",
    "n_markers_used",
    "n_references",
    "n_samples_polymorphic_data",
    "mean_pct_polymorphic",
    "max_pct_polymorphic",
    "sample_max_pct_polymorphic",
]


def main():
    parser = argparse.ArgumentParser(
        description="Summarize StrainPhlAn4 output across all clades.")
    parser.add_argument("--root", required=True,
                        help="Root dir containing t__SGB* subdirectories")
    parser.add_argument("--names", default=None,
                        help="Path to output_sgb_names.tsv (optional but recommended)")
    parser.add_argument("--out", default="strainphlan_summary.tsv",
                        help="Output TSV (default: strainphlan_summary.tsv)")
    parser.add_argument("--debug-info", metavar="CLADE",
                        help="Print .info file for one clade and exit")
    args = parser.parse_args()

    root = Path(args.root)

    # Debug mode
    if args.debug_info:
        info_path = root / args.debug_info / f"{args.debug_info}.info"
        if info_path.exists():
            lines = info_path.read_text(errors="replace").splitlines()
            print(f"=== {info_path} (first 60 lines) ===")
            for line in lines[:60]:
                print(line)
        else:
            print(f"Not found: {info_path}")
        return

    if not root.is_dir():
        raise SystemExit(f"ERROR: --root '{root}' is not a directory.")

    # Load names table
    sgb_names = load_sgb_names(args.names)
    if sgb_names:
        print(f"Loaded {len(sgb_names)} SGB names from {args.names}")
    else:
        print("No SGB names file loaded (species/genus columns will be empty).")

    clade_dirs = sorted([d for d in root.iterdir()
                         if d.is_dir() and d.name.startswith("t__SGB")])
    if not clade_dirs:
        raise SystemExit(f"ERROR: No t__SGB* subdirectories found in '{root}'.")

    print(f"Found {len(clade_dirs)} clade directories. Summarizing...")

    rows = []
    for i, clade_dir in enumerate(clade_dirs, 1):
        try:
            row = summarize_clade(clade_dir, sgb_names)
        except Exception as e:
            print(f"  WARNING: Failed to parse {clade_dir.name}: {e}")
            row = {f: None for f in FIELDNAMES}
            row["clade"] = clade_dir.name
        rows.append(row)
        if i % 100 == 0:
            print(f"  Processed {i}/{len(clade_dirs)}...")

    out_path = Path(args.out)
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nDone! Written to: {out_path}")
    parsed_tree    = sum(1 for r in rows if r["n_samples_in_tree"] is not None)
    parsed_markers = sum(1 for r in rows if r["n_markers_used"] is not None)
    parsed_poly    = sum(1 for r in rows if r["n_samples_polymorphic_data"] > 0)
    named          = sum(1 for r in rows if r["species"])
    print(f"  Tree tips parsed:          {parsed_tree}/{len(rows)}")
    print(f"  Markers from .info parsed: {parsed_markers}/{len(rows)}")
    print(f"  Polymorphic data parsed:   {parsed_poly}/{len(rows)}")
    print(f"  Species names matched:     {named}/{len(rows)}")


if __name__ == "__main__":
    main()