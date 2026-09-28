#!/usr/bin/env python3
"""Generate nf-core/chipseq (including DAP-seq) and nf-core/atacseq samplesheets.

Accept a header table (sample_id, antibody, control_id, replicate, fastq_path),
legacy positional columns, or explicit fastq_1 and fastq_2 columns.
File-level comments can specify #@ pair_pattern: _R{R}_001.fastq.gz
and #@ single_end: true. See `chromaflow help samples` for examples.
"""

import argparse
import csv
import glob
import os
import re
import sys
from collections import defaultdict

# ── FASTQ 扩展名（按长优先匹配）──
FASTQ_EXTS = (
    ".fastq.gz", ".fq.gz", ".fastq.bz2", ".fq.bz2", ".fastq.xz", ".fq.xz",
    ".fastq.zst", ".fq.zst", ".fastq", ".fq",
)

# ── 双端记号：R1/R2、r1/r2、read1/read2、1/2，可带 P/U 后缀(trimmomatic) ──
MATE_TOKEN = re.compile(r"^(?P<tag>R|r|read|Read|READ)?(?P<num>[12])(?P<suf>[PU])?$", re.ASCII)
# 无分隔符结尾的情况: sampleR1.fq.gz
MATE_TAIL = re.compile(r"(?P<tag>[Rr]|read|Read|READ)(?P<num>[12])$", re.ASCII)

SPLIT_TOKENS = re.compile(r"([._-])")


# ────────────────────────────── FASTQ 配对 ──────────────────────────────
def split_fastq_ext(name):
    """Split a FASTQ extension; return (name, '') for other file types."""
    low = name.lower()
    for ext in FASTQ_EXTS:
        if low.endswith(ext):
            return name[: -len(ext)], name[-len(ext):]
    return name, ""


def mate_of(path, pattern=None):
    """Infer the R2 path from an R1 file.

    Use a user template such as '_good_{R}.fq.gz', or detect the rightmost mate token.
    Return (r2_path, token); return (None, None) if no mate can be inferred.
    """
    d, base = os.path.split(path)

    # ── 1) 用户显式模板 ──
    if pattern:
        r1_suf, r2_suf = pattern.replace("{R}", "1"), pattern.replace("{R}", "2")
        if base.endswith(r1_suf):
            return os.path.join(d, base[: -len(r1_suf)] + r2_suf), pattern
        return None, None

    # ── 2) 自动：按分隔符切 token，从右往左找第一个双端记号 ──
    stem, ext = split_fastq_ext(base)
    if not ext:
        return None, None
    parts = SPLIT_TOKENS.split(stem)          # ['X', '_', 'R1', '_', '001']
    for i in range(len(parts) - 1, -1, -1):
        if i % 2 == 1:                        # 分隔符位
            continue
        m = MATE_TOKEN.match(parts[i])
        if m and m.group("num") == "1":
            np = list(parts)
            np[i] = "%s2%s" % (m.group("tag") or "", m.group("suf") or "")
            return os.path.join(d, "".join(np) + ext), parts[i]
        if m and m.group("num") == "2":
            return None, None                 # 这是 R2，跳过（由 R1 侧Pair）

    # ── 3) 无分隔符收尾: sampleR1.fq.gz ──
    m = MATE_TAIL.search(stem)
    if m and m.group("num") == "1":
        return os.path.join(d, stem[: m.start()] + m.group("tag") + "2" + ext), m.group(0)
    return None, None


def is_mate2(path, pattern=None):
    """Identify R2 files so they are not treated as independent single-end samples."""
    d, base = os.path.split(path)
    if pattern:
        return base.endswith(pattern.replace("{R}", "2"))
    stem, ext = split_fastq_ext(base)
    if not ext:
        return False
    parts = SPLIT_TOKENS.split(stem)
    for i in range(len(parts) - 1, -1, -1):
        if i % 2 == 1:
            continue
        m = MATE_TOKEN.match(parts[i])
        if m:
            return m.group("num") == "2"
    m = MATE_TAIL.search(stem)
    return bool(m and m.group("num") == "2")


def collect_candidates(path):
    """Accept a directory, filename prefix, glob, or a single FASTQ file."""
    path = path.strip()
    cands = []
    if os.path.isfile(path):
        cands = [path]
    elif os.path.isdir(path):
        for ext in FASTQ_EXTS:
            cands += glob.glob(os.path.join(path, "*" + ext))
            cands += glob.glob(os.path.join(path, "**", "*" + ext), recursive=True)
    elif any(c in path for c in "*?["):
        cands = [f for f in glob.glob(path, recursive=True) if split_fastq_ext(f)[1]]
    else:
        parent = os.path.dirname(path) or "."
        prefix = os.path.basename(path)
        if os.path.isdir(parent):
            cands = [os.path.join(parent, f) for f in os.listdir(parent)
                     if f.startswith(prefix) and split_fastq_ext(f)[1]]
    return sorted(set(os.path.abspath(f) for f in cands))


def find_fastqs(path, pattern=None, single_end=False, verbose=True):
    """Find FASTQ pairs; use an empty fastq_2 for single-end reads."""
    cands = collect_candidates(path)
    if not cands:
        print("  [WARN] FASTQ not found: %s" % path, file=sys.stderr)
        return []

    if single_end:
        return [(f, "") for f in cands]

    pairs, used, tokens = [], set(), set()
    for fq in cands:
        if fq in used:
            continue
        r2, token = mate_of(fq, pattern)
        if r2 and os.path.exists(r2):
            pairs.append((fq, r2))
            used.update([fq, r2])
            if token:
                tokens.add(token)
        elif r2:
            print("  [WARN] Mate file missing; treating as single-end: %s (expected %s)"
                  % (os.path.basename(fq), os.path.basename(r2)), file=sys.stderr)
            pairs.append((fq, ""))
            used.add(fq)
        elif not is_mate2(fq, pattern):
            pairs.append((fq, ""))
            used.add(fq)
    if verbose and tokens:
        print("    · Pairing tokens: %s" % ", ".join(sorted(tokens)))
    return pairs


# ────────────────────────────── 样本记录解析 ──────────────────────────────
COLUMN_ALIASES = {
    "sample_id": "sample_id", "sample": "sample_id", "sampleid": "sample_id", "name": "sample_id",
    "antibody": "antibody", "ab": "antibody", "target": "antibody",
    "control_id": "control_id", "control": "control_id", "input": "control_id",
    "control_sample": "control_id",
    "replicate": "replicate", "rep": "replicate", "replicate_id": "replicate",
    "fastq_path": "fastq_path", "fastq": "fastq_path", "path": "fastq_path",
    "dir": "fastq_path", "fastq_dir": "fastq_path", "reads": "fastq_path",
    "fastq_1": "fastq_1", "fq1": "fastq_1", "r1": "fastq_1", "read1": "fastq_1",
    "fastq_2": "fastq_2", "fq2": "fastq_2", "r2": "fastq_2", "read2": "fastq_2",
    "pair_pattern": "pair_pattern", "pattern": "pair_pattern",
    "single_end": "single_end", "se": "single_end",
}
POSITIONAL = ["sample_id", "antibody", "control_id", "replicate", "fastq_path"]
EMPTY = ("", "-", "NA", "na", ".", "none")
DIRECTIVE = re.compile(r"^#\s*@?\s*([A-Za-z_]+)\s*[:=]\s*(.+?)\s*$")


def _split(line):
    if "\t" in line:
        return [c.strip() for c in line.split("\t")]
    if "," in line and line.count(",") >= line.count(" "):
        return [c.strip() for c in line.split(",")]
    return line.split()


def _is_header(cols):
    known = sum(1 for c in cols if c.strip().lower().lstrip("#") in COLUMN_ALIASES)
    return known >= max(2, len(cols) - 1)


def _truthy(v):
    return str(v).strip().lower() in ("true", "yes", "1", "y")


def parse_samples(filepath):
    """Return a list of record dictionaries and file-level directives."""
    records, directives, header = [], {}, None
    with open(filepath, encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n").strip()
            if not line:
                continue
            if line.startswith("#"):
                m = DIRECTIVE.match(line)
                if m and m.group(1).lower() in ("pair_pattern", "pattern", "single_end", "se"):
                    directives[COLUMN_ALIASES[m.group(1).lower()]] = m.group(2)
                elif m and _is_header(_split(line.lstrip("#").strip())) and header is None:
                    header = [COLUMN_ALIASES.get(c.lower(), c.lower())
                              for c in _split(line.lstrip("#").strip())]
                continue

            cols = _split(line)
            if header is None and _is_header(cols):
                header = [COLUMN_ALIASES.get(c.lower(), c.lower()) for c in cols]
                continue

            if header:
                rec = {header[i]: cols[i] for i in range(min(len(header), len(cols)))}
            else:
                if len(cols) < 4:
                    print("  [ERROR] Line %d has too few columns (positional format needs >=4, or use a header): %s"
                          % (lineno, line), file=sys.stderr)
                    continue
                rec = {POSITIONAL[i]: cols[i] for i in range(min(len(POSITIONAL), len(cols)))}

            if not rec.get("sample_id"):
                print("  [ERROR] Line %d is missing sample_id" % lineno, file=sys.stderr)
                continue
            for k in ("antibody", "control_id", "fastq_path", "fastq_1", "fastq_2"):
                if rec.get(k, "").strip() in EMPTY:
                    rec[k] = ""
            rec.setdefault("replicate", "1")
            if not str(rec["replicate"]).strip() or rec["replicate"] in EMPTY:
                rec["replicate"] = "1"
            rec["_line"] = lineno
            records.append(rec)
    return records, directives


def resolve_reads(rec, pattern, single_end):
    """Expand one record into FASTQ pairs."""
    if rec.get("fastq_1"):
        return [(os.path.abspath(rec["fastq_1"]),
                 os.path.abspath(rec["fastq_2"]) if rec.get("fastq_2") else "")]
    if not rec.get("fastq_path"):
        return [("", "")]
    pat = rec.get("pair_pattern") or pattern
    se = _truthy(rec["single_end"]) if rec.get("single_end") else single_end
    return find_fastqs(rec["fastq_path"], pat, se)


# ────────────────────────────── samplesheet 生成 ──────────────────────────────
def build_rows(records, pipeline, pattern, single_end):
    ctrl_reps = defaultdict(list)
    for r in records:
        if not r["antibody"]:
            ctrl_reps[r["sample_id"]].append(r["replicate"])

    rows = []
    for r in records:
        reads = resolve_reads(r, pattern, single_end)
        n_pe = sum(1 for _, f2 in reads if f2)
        print("  · %-18s rep%-3s %d libraries (%d paired-end / %d single-end)"
              % (r["sample_id"], r["replicate"], len(reads), n_pe, len(reads) - n_pe))
        for f1, f2 in reads:
            row = {"sample": r["sample_id"], "fastq_1": f1, "fastq_2": f2,
                   "replicate": r["replicate"]}
            if pipeline == "chipseq":
                cr = ""
                if r["control_id"] and r["antibody"]:
                    avail = ctrl_reps.get(r["control_id"], [])
                    cr = r["replicate"] if r["replicate"] in avail else (avail[0] if avail else "")
                row.update({"antibody": r["antibody"], "control": r["control_id"],
                            "control_replicate": cr})
            rows.append(row)
    return rows


def write_csv(path, pipeline, rows):
    header = (["sample", "fastq_1", "fastq_2", "replicate", "antibody", "control", "control_replicate"]
              if pipeline == "chipseq" else ["sample", "fastq_1", "fastq_2", "replicate"])
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def validate(path, pipeline):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    issues, sr = [], defaultdict(set)
    for row in rows:
        try:
            sr[row["sample"]].add(int(row["replicate"]))
        except (ValueError, TypeError):
            issues.append("  ⚠  '%s' replicate is not an integer: %r" % (row["sample"], row["replicate"]))
        for k in ("fastq_1", "fastq_2"):
            if row.get(k) and not os.path.exists(row[k]):
                issues.append("  ⚠  File does not exist: %s" % row[k])
        if not row.get("fastq_1"):
            issues.append("  ⚠  '%s' has no fastq_1; check fastq_path" % row["sample"])
    for s, reps in sr.items():
        if 1 not in reps:
            issues.append("  ⚠  '%s' replicates do not include 1 (nf-core requires numbering from 1)" % s)
    if pipeline == "chipseq":
        all_s = set(sr.keys())
        for row in rows:
            c = row.get("control", "")
            if c and c not in all_s:
                issues.append("  ⚠  '%s' references missing control '%s'" % (row["sample"], c))
    if issues:
        print("\n⚠  %d issues:" % len(issues))
        for i in dict.fromkeys(issues):
            print(i)
    else:
        print("\n✅ Validation passed")
    return not issues


def preview(path):
    with open(path) as f:
        lines = f.read().strip().split("\n")
    print("\n%s\n📋 Preview (%d rows):\n%s" % ("─" * 90, len(lines) - 1, "─" * 90))
    print(lines[0])
    for l in lines[1:]:
        parts = l.split(",")
        for i in (1, 2):
            if i < len(parts) and len(parts[i]) > 50:
                parts[i] = "..." + parts[i][-47:]
        print(",".join(parts))
    print("─" * 90)


def main():
    p = argparse.ArgumentParser(description="Generate an nf-core sample sheet")
    p.add_argument("samples_txt")
    p.add_argument("-o", "--output", default="samplesheet.csv")
    p.add_argument("-p", "--pipeline", choices=["chipseq", "atacseq"], default="chipseq")
    p.add_argument("--pair-pattern", default=None,
                   help="Read-pair naming pattern with {R}=1/2, e.g. '_good_{R}.fq.gz' (default: automatic detection)")
    p.add_argument("--single-end", action="store_true", help="Treat all inputs as single-end")
    p.add_argument("--list-fastq", action="store_true", help="Print paired FASTQ files without writing a samplesheet")
    p.add_argument("--strict", action="store_true", help="Exit nonzero if validation finds issues")
    p.add_argument("--validate-only", action="store_true")
    args = p.parse_args()

    if args.validate_only:
        sys.exit(0 if validate(args.samples_txt, args.pipeline) or not args.strict else 1)

    print("📄 Reading: %s" % args.samples_txt)
    recs, directives = parse_samples(args.samples_txt)
    if not recs:
        print("[ERROR] No valid records", file=sys.stderr)
        sys.exit(1)

    pattern = args.pair_pattern or directives.get("pair_pattern") or None
    single_end = args.single_end or _truthy(directives.get("single_end", "false"))
    print("   %d records | Pairing mode: %s%s"
          % (len(recs), ("pattern " + pattern) if pattern else "automatic detection",
             " | forced single-end" if single_end else ""))
    print("🔧 Generating nf-core/%s sample sheet..." % args.pipeline)

    rows = build_rows(recs, args.pipeline, pattern, single_end)

    if args.list_fastq:
        for r in rows:
            print("%s\trep%s\t%s\t%s" % (r["sample"], r["replicate"], r["fastq_1"], r["fastq_2"]))
        return

    write_csv(args.output, args.pipeline, rows)
    print("\n✅ Generated: %s (%d rows)" % (args.output, len(rows)))
    preview(args.output)
    ok = validate(args.output, args.pipeline)
    if args.strict and not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
