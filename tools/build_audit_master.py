#!/usr/bin/env python3
"""
build_audit_master.py — regenerate metadata/audit_master.tsv from tools/tier_audit.py.

For every bundle (AT, CR, MP) it runs the tier audit against the upstream folder and
writes one row per .m file with:

  file_name  bundle_id  tier  provenance  invoked  module_block  role_1line  rel_path

`tier` follows docs/tiers.md (T4 = present but not invoked, precedence over T1-T3).
`module_block` and `role_1line` are carried over from the previous audit_master.tsv
when the file name is known; otherwise the block is inferred from the name and the role
is the first comment line of the file. Rows for upstream files that were never ported
(NOT_PORTED) are kept as they were.

Usage
    python tools/build_audit_master.py --upstream /path/to/EpiProfile2.0_1Basic
"""
import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MASTER = REPO / "metadata" / "audit_master.tsv"
BUNDLES = ["AT", "CR", "MP"]

BLOCK_RULES = [
    (r"^(H3|H4|HH1|HH2A|HH2B|H2A|H2B)[A-Za-z]*_", "histone_module"),
    (r"^(EpiProfile|DrawISOProfile|run_calibrate)", "runner"),
    (r"^(Output|output_)", "io_output"),
    (r"^(ReadInput|check_otherparas|Raw2MS)", "io_control"),
    (r"^(GetMS1ScanNo|GetMS2ScanNo)", "io_input"),
    (r"^(get_histone|get_rts|get_area|get_theo|get_key|get_main|get_mod|find_pair|find_triple|GetProfiles|GetTopBottom|GetIPV|GetLocal|JudgeLocalmaxmin|relocateD)", "histone_quant"),
    (r"^(calculate_pepmz|Getaamass|GetMods|check_ref|check_layout|GetPSM|GetBenchmark)", "core_quant"),
    (r"^(draw_layout|OutputFigures)", "visualization"),
]


def infer_block(name: str) -> str:
    for pat, block in BLOCK_RULES:
        if re.match(pat, name):
            return block
    return "other"


def first_comment(path: Path) -> str:
    try:
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines()[:12]:
            s = line.strip()
            if s.startswith("%"):
                s = s.lstrip("%").strip()
                if s and s != "%":
                    return re.sub(r"\s+", " ", s)[:120]
    except OSError:
        pass
    return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--upstream", required=True, type=Path)
    ap.add_argument("--out", type=Path, default=MASTER)
    args = ap.parse_args()

    old = {}
    not_ported = []
    if MASTER.exists():
        with open(MASTER, encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("tier") == "NOT_PORTED":
                    not_ported.append(row)
                else:
                    old.setdefault(row["file_name"], row)

    rows = []
    for b in BUNDLES:
        src = REPO / "bundles" / b / "src"
        if not src.is_dir():
            continue
        tmp = REPO / "bundles" / b / "metadata" / f"tier_audit_{b}.tsv"
        tmp.parent.mkdir(exist_ok=True)
        subprocess.run([sys.executable, str(REPO / "tools" / "tier_audit.py"), "--bundle", str(src),
                        "--upstream", str(args.upstream), "--bundle-id", b, "--out", str(tmp)],
                       check=True)
        with open(tmp, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                o = old.get(r["file_name"], {})
                role = o.get("role_1line") or first_comment(src / r["rel_path"])
                block = o.get("module_block") or infer_block(r["file_name"])
                if r["tier"] == "T4" and "not invoked" not in role.lower():
                    role = (role + " " if role else "") + "[not invoked by the bundle]"
                rows.append({"file_name": r["file_name"], "bundle_id": b, "tier": r["tier"],
                             "provenance": r["provenance"], "invoked": r["invoked"],
                             "module_block": block, "role_1line": role, "rel_path": r["rel_path"]})
    for r in not_ported:
        rows.append({"file_name": r["file_name"], "bundle_id": r.get("bundle_id", "AT"), "tier": "NOT_PORTED",
                     "provenance": "upstream", "invoked": "no", "module_block": r.get("module_block", ""),
                     "role_1line": r.get("role_1line", ""), "rel_path": ""})

    cols = ["file_name", "bundle_id", "tier", "provenance", "invoked", "module_block", "role_1line", "rel_path"]
    text = "\t".join(cols) + "\n" + "\n".join("\t".join(r[c] for c in cols) for r in rows) + "\n"
    args.out.write_text(text, encoding="utf-8", newline="\n")
    summary = {}
    for r in rows:
        summary.setdefault(r["bundle_id"], {}).setdefault(r["tier"], 0)
        summary[r["bundle_id"]][r["tier"]] += 1
    print(f"wrote {args.out} ({len(rows)} rows)")
    for b, c in summary.items():
        print(f"  {b}: " + ", ".join(f"{k}={v}" for k, v in sorted(c.items())))


if __name__ == "__main__":
    main()
