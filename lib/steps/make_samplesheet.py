#!/usr/bin/env python3
"""
make_samplesheet.py — 从 samples 记录表生成 nf-core samplesheet

支持: nf-core/chipseq v2.1.0 (含 DAP-seq) | nf-core/atacseq v2.1.2

样本记录表（两种写法，可混用注释）:

  1) 带表头（推荐，列可任意顺序、可缺省、可加自定义列）:
       sample_id  antibody  control_id  replicate  fastq_path
       Dof_ChIP   crDOF     Input       1          /data/raw/Dof_rep1

  2) 传统位置列（向后兼容）:
       sample_id <TAB> antibody <TAB> control_id <TAB> replicate <TAB> fastq_path

  也可直接给出双端文件（跳过自动配对）:
       sample_id  replicate  fastq_1              fastq_2

指令行（写在文件顶部的注释里，对整份表生效）:
       #@ pair_pattern: _good_{R}.fq.gz     # {R} 处填 1/2
       #@ single_end: true                  # 全部按单端处理

用法:
  python make_samplesheet.py samples.txt -o samplesheet.csv --pipeline chipseq
  python make_samplesheet.py samples.txt -o out.csv --pair-pattern '_R{R}_001.fastq.gz'
  python make_samplesheet.py samples.txt --list-fastq        # 只看配对结果，不写文件
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
    """('X_R1_001.fastq.gz') → ('X_R1_001', '.fastq.gz')；非 fastq 返回 (name, '')"""
    low = name.lower()
    for ext in FASTQ_EXTS:
        if low.endswith(ext):
            return name[: -len(ext)], name[-len(ext):]
    return name, ""


def mate_of(path, pattern=None):
    """给定 R1 文件路径，推断 R2 路径。

    pattern: 用户模板，如 '_good_{R}.fq.gz'；None 则自动识别最右侧的双端记号。
    返回 (r2_path, 记号描述)；无法判定返回 (None, None)。
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
            return None, None                 # 这是 R2，跳过（由 R1 侧配对）

    # ── 3) 无分隔符收尾: sampleR1.fq.gz ──
    m = MATE_TAIL.search(stem)
    if m and m.group("num") == "1":
        return os.path.join(d, stem[: m.start()] + m.group("tag") + "2" + ext), m.group(0)
    return None, None


def is_mate2(path, pattern=None):
    """判断某文件是否是 R2（避免把 R2 当成独立单端样本）。"""
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
    """path 可为目录 / 文件前缀 / glob / 单个 fastq 文件。"""
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
    """发现并配对 FASTQ；返回 [(fastq_1, fastq_2), ...]，单端时 fastq_2 为空。"""
    cands = collect_candidates(path)
    if not cands:
        print("  [WARN] 未找到 fastq: %s" % path, file=sys.stderr)
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
            print("  [WARN] 缺配对文件，按单端处理: %s (期望 %s)"
                  % (os.path.basename(fq), os.path.basename(r2)), file=sys.stderr)
            pairs.append((fq, ""))
            used.add(fq)
        elif not is_mate2(fq, pattern):
            pairs.append((fq, ""))
            used.add(fq)
    if verbose and tokens:
        print("    · 配对记号: %s" % ", ".join(sorted(tokens)))
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
    """返回 (records, directives)。records 为 dict 列表。"""
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
                    print("  [ERROR] 行%d 列数不足（位置格式需 ≥4 列，或改用表头格式）: %s"
                          % (lineno, line), file=sys.stderr)
                    continue
                rec = {POSITIONAL[i]: cols[i] for i in range(min(len(POSITIONAL), len(cols)))}

            if not rec.get("sample_id"):
                print("  [ERROR] 行%d 缺 sample_id" % lineno, file=sys.stderr)
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
    """一条记录 → [(fastq_1, fastq_2), ...]"""
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
        print("  · %-18s rep%-3s %d 个文件 (%d 双端 / %d 单端)"
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
            issues.append("  ⚠  '%s' replicate 非整数: %r" % (row["sample"], row["replicate"]))
        for k in ("fastq_1", "fastq_2"):
            if row.get(k) and not os.path.exists(row[k]):
                issues.append("  ⚠  文件不存在: %s" % row[k])
        if not row.get("fastq_1"):
            issues.append("  ⚠  '%s' 无 fastq_1（fastq_path 没找到文件？）" % row["sample"])
    for s, reps in sr.items():
        if 1 not in reps:
            issues.append("  ⚠  '%s' replicate 不含 1（nf-core 要求从 1 开始）" % s)
    if pipeline == "chipseq":
        all_s = set(sr.keys())
        for row in rows:
            c = row.get("control", "")
            if c and c not in all_s:
                issues.append("  ⚠  '%s' 引用了不存在的 control '%s'" % (row["sample"], c))
    if issues:
        print("\n⚠  %d 个问题:" % len(issues))
        for i in dict.fromkeys(issues):
            print(i)
    else:
        print("\n✅ 验证通过")
    return not issues


def preview(path):
    with open(path) as f:
        lines = f.read().strip().split("\n")
    print("\n%s\n📋 预览 (%d 行):\n%s" % ("─" * 90, len(lines) - 1, "─" * 90))
    print(lines[0])
    for l in lines[1:]:
        parts = l.split(",")
        for i in (1, 2):
            if i < len(parts) and len(parts[i]) > 50:
                parts[i] = "..." + parts[i][-47:]
        print(",".join(parts))
    print("─" * 90)


def main():
    p = argparse.ArgumentParser(description="生成 nf-core samplesheet")
    p.add_argument("samples_txt")
    p.add_argument("-o", "--output", default="samplesheet.csv")
    p.add_argument("-p", "--pipeline", choices=["chipseq", "atacseq"], default="chipseq")
    p.add_argument("--pair-pattern", default=None,
                   help="双端命名模板，{R} 处为 1/2，如 '_good_{R}.fq.gz'（默认自动识别）")
    p.add_argument("--single-end", action="store_true", help="全部按单端处理")
    p.add_argument("--list-fastq", action="store_true", help="只打印配对结果，不写 samplesheet")
    p.add_argument("--strict", action="store_true", help="有问题时以非 0 退出")
    p.add_argument("--validate-only", action="store_true")
    args = p.parse_args()

    if args.validate_only:
        sys.exit(0 if validate(args.samples_txt, args.pipeline) or not args.strict else 1)

    print("📄 读取: %s" % args.samples_txt)
    recs, directives = parse_samples(args.samples_txt)
    if not recs:
        print("[ERROR] 无有效记录", file=sys.stderr)
        sys.exit(1)

    pattern = args.pair_pattern or directives.get("pair_pattern") or None
    single_end = args.single_end or _truthy(directives.get("single_end", "false"))
    print("   %d 条记录 | 配对方式: %s%s"
          % (len(recs), ("模板 " + pattern) if pattern else "自动识别",
             " | 强制单端" if single_end else ""))
    print("🔧 生成 nf-core/%s samplesheet..." % args.pipeline)

    rows = build_rows(recs, args.pipeline, pattern, single_end)

    if args.list_fastq:
        for r in rows:
            print("%s\trep%s\t%s\t%s" % (r["sample"], r["replicate"], r["fastq_1"], r["fastq_2"]))
        return

    write_csv(args.output, args.pipeline, rows)
    print("\n✅ 已生成: %s (%d 行)" % (args.output, len(rows)))
    preview(args.output)
    ok = validate(args.output, args.pipeline)
    if args.strict and not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
