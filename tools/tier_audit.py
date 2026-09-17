#!/usr/bin/env python3
"""
tier_audit.py — provenance (T1/T2/T3) and invocation (T4) audit of a bundle.

Implements the rules of docs/tiers.md:

  T1  file identical to the upstream file of the same name (content normalised
      for line endings and trailing whitespace)
  T2  same name exists upstream, content differs
  T3  name does not exist upstream (new in the PLANTS family)
  T4  present in the bundle but never invoked by another bundle file
      (T4 takes precedence over T1-T3, so T1+T2+T3+T4 = number of .m files)

Invocation = the file's basename appears as a word in some OTHER .m file of the
bundle, ignoring comments (`% ...`) and string literals. Entry points named
with --entry (default: EpiProfile, run_calibrate) count as invoked.

Usage
    python tools/tier_audit.py --bundle bundles/AT/src \
        --upstream /path/to/EpiProfile2.0_1Basic --bundle-id AT \
        --out bundles/AT/metadata/tier_audit_AT.tsv

The upstream folder is not part of this repository (EpiProfile 2.0 basic,
Yuan et al. 2018); point --upstream at your copy.
"""
import argparse
import hashlib
import re
import sys
from pathlib import Path


def norm(text: str) -> str:
    lines = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    lines = [ln.rstrip() for ln in lines]
    while lines and lines[-1] == "":
        lines.pop()
    return "\n".join(lines) + "\n"


def read(p: Path) -> str:
    return p.read_text(encoding="utf-8", errors="replace")


def strip_comments_and_strings(code: str) -> str:
    """Remove %-comments and '...' string literals (keeps transposes)."""
    out = []
    for line in code.split("\n"):
        res = []
        in_str = False
        i = 0
        prev = ""
        while i < len(line):
            ch = line[i]
            if in_str:
                if ch == "'":
                    if i + 1 < len(line) and line[i + 1] == "'":   # escaped ''
                        i += 2
                        continue
                    in_str = False
                i += 1
                continue
            if ch == "%":
                break
            if ch == "'":
                # transpose if preceded by identifier char, ')' , ']', '}', '.', or "'"
                if prev and (prev.isalnum() or prev in ")]}._'"):
                    res.append(ch)
                else:
                    in_str = True
                i += 1
                prev = ch
                continue
            res.append(ch)
            prev = ch if not ch.isspace() else prev
            i += 1
        out.append("".join(res))
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", required=True, type=Path, help="bundle src folder (searched recursively)")
    ap.add_argument("--upstream", type=Path, default=None, help="upstream EpiProfile 2.0 basic folder")
    ap.add_argument("--bundle-id", default="AT")
    ap.add_argument("--entry", nargs="*", default=["EpiProfile", "run_calibrate"],
                    help="entry points that count as invoked")
    ap.add_argument("--out", type=Path, default=None, help="TSV output (default: stdout)")
    args = ap.parse_args()

    files = sorted(args.bundle.rglob("*.m"))
    if not files:
        sys.exit(f"no .m files under {args.bundle}")
    names = [f.stem for f in files]
    dup = {n for n in names if names.count(n) > 1}

    upstream = {}
    if args.upstream:
        for f in sorted(args.upstream.rglob("*.m")):
            upstream.setdefault(f.stem, f)

    # code bodies without comments/strings, for invocation search
    bodies = {f: strip_comments_and_strings(read(f)) for f in files}
    patterns = {n: re.compile(r"(?<![A-Za-z0-9_])" + re.escape(n) + r"(?![A-Za-z0-9_])") for n in set(names)}

    rows = []
    counts = {"T1": 0, "T2": 0, "T3": 0, "T4": 0}
    for f in files:
        n = f.stem
        h = hashlib.sha256(norm(read(f)).encode()).hexdigest()[:12]
        if n in upstream:
            prov = "T1" if norm(read(upstream[n])) == norm(read(f)) else "T2"
            # relative to --upstream, so the TSV carries no workstation path
            up = str(upstream[n].relative_to(args.upstream)).replace("\\", "/")
        else:
            prov = "T3"
            up = ""
        callers = [g.stem for g in files if g != f and patterns[n].search(bodies[g])]
        invoked = bool(callers) or n in args.entry
        tier = prov if invoked else "T4"
        counts[tier] += 1
        notes = []
        if n in dup:
            notes.append("DUPLICATE name in bundle (shadowing)")
        if n in args.entry:
            notes.append("entry point")
        rows.append({
            "file_name": f.name, "bundle_id": args.bundle_id, "tier": tier, "provenance": prov,
            "invoked": "yes" if invoked else "no", "invoked_by": ",".join(sorted(set(callers))),
            "rel_path": str(f.relative_to(args.bundle)).replace("\\", "/"),
            "sha256_12": h, "upstream_file": up, "notes": "; ".join(notes),
        })

    cols = ["file_name", "bundle_id", "tier", "provenance", "invoked", "invoked_by", "rel_path",
            "sha256_12", "upstream_file", "notes"]
    lines = ["\t".join(cols)] + ["\t".join(r[c] for c in cols) for r in rows]
    text = "\n".join(lines) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text, encoding="utf-8", newline="\n")
    else:
        sys.stdout.write(text)
    total = len(files)
    print(f"# {args.bundle_id}: {total} .m  T1={counts['T1']} T2={counts['T2']} T3={counts['T3']} "
          f"T4={counts['T4']}  (sum={sum(counts.values())})"
          + (f"  upstream={len(upstream)} files" if upstream else "  (no upstream given: T1/T2 not resolved)"),
          file=sys.stderr)
    if dup:
        print(f"# duplicate names: {sorted(dup)}", file=sys.stderr)


if __name__ == "__main__":
    main()
