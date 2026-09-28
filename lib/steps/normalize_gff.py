#!/usr/bin/env python3
"""Conservative offline GFF3 normalization. Existing records are immutable.

Unknown ontology terms are preserved, never mistaken for missing graph nodes.
The bundled role table is finite; local extensions are explicit.
"""
import argparse
import bz2
from collections import Counter, defaultdict, deque
from contextlib import contextmanager
import copy
import gzip
import hashlib
import io
import json
import lzma
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import zipfile
from urllib.parse import quote, unquote

VERSION = "1.2.3"
POLICY = "conservative-v1"
ROLE_FILE = Path(__file__).with_name("gff_roles.json")


class GFFError(ValueError):
    pass


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


@contextmanager
def opener(path):
    with open(path, "rb") as handle:
        magic = handle.read(6)
    factory = (gzip.open if magic.startswith(b"\x1f\x8b") else
               bz2.open if magic.startswith(b"BZh") else
               lzma.open if magic.startswith(b"\xfd7zXZ\x00") else None)
    if factory:
        with factory(path, "rt", encoding="utf-8") as handle:
            yield handle
    elif magic.startswith(b"PK\x03\x04"):
        with zipfile.ZipFile(path) as archive:
            members = [m for m in archive.infolist() if not m.is_dir()]
            if len(members) != 1:
                raise GFFError("GFF zip input must contain exactly one file")
            with archive.open(members[0]) as raw, io.TextIOWrapper(raw, encoding="utf-8") as handle:
                yield handle
    elif magic.startswith(b"\x28\xb5\x2f\xfd"):
        with tempfile.TemporaryFile() as errors:
            proc = subprocess.Popen(["zstd", "-dc", "--", str(path)], stdout=subprocess.PIPE, stderr=errors)
            try:
                with io.TextIOWrapper(proc.stdout, encoding="utf-8") as handle:
                    yield handle
                if proc.wait():
                    raise GFFError("zstd decompression failed")
            finally:
                if proc.poll() is None:
                    proc.terminate()
                    proc.wait()
    else:
        with open(path, encoding="utf-8") as handle:
            yield handle


def decode(value):
    if re.search(r"%(?![0-9a-fA-F]{2})", value):
        raise GFFError("invalid percent escape in attribute")
    return unquote(value, encoding="utf-8", errors="strict")


def encode(value):
    return quote(value, safe=".:^*$@!+_?-|")


class Feature:
    def __init__(self, cols, line_number=0):
        if len(cols) != 9:
            raise GFFError("expected exactly 9 tab-separated columns")
        self.cols = tuple(cols)
        self.line_number = line_number
        self.seqid, self.source, self.type = cols[:3]
        try:
            self.start, self.end = int(cols[3]), int(cols[4])
        except ValueError as exc:
            raise GFFError("coordinates must be integers") from exc
        if self.start < 1 or self.end < self.start:
            raise GFFError("expected 1 <= start <= end")
        self.strand, self.phase = cols[6:8]
        if self.strand not in {"+", "-", ".", "?"} or self.phase not in {"0", "1", "2", "."}:
            raise GFFError("invalid strand or phase")
        if cols[5] != ".":
            try:
                if not math.isfinite(float(cols[5])):
                    raise ValueError()
            except ValueError as exc:
                raise GFFError("invalid score") from exc
        if not self.seqid or self.seqid == "." or not self.type or self.type == ".":
            raise GFFError("seqid and type are required")
        self.attrs = {}
        if cols[8] != ".":
            for item in cols[8].rstrip(";").split(";"):
                if "=" not in item:
                    raise GFFError("expected GFF3 key=value attributes (not GTF)")
                key, value = item.split("=", 1)
                if not key or key in self.attrs:
                    raise GFFError("empty or duplicate attribute key: " + key)
                decode(value)
                self.attrs[key] = value
        raw_id = self.attrs.get("ID", "")
        if "ID" in self.attrs and (not raw_id or "," in raw_id):
            raise GFFError("ID must contain one nonempty escaped identifier")
        self.id = decode(raw_id)
        raw_parent = self.attrs.get("Parent")
        if raw_parent is not None and any(not p for p in raw_parent.split(",")):
            raise GFFError("Parent contains an empty identifier")
        # Split BEFORE decoding: %2C is a literal comma in a single ID.
        self.parents = tuple(decode(p) for p in raw_parent.split(",")) if raw_parent else ()
        if len(set(self.parents)) != len(self.parents):
            raise GFFError("duplicate Parent edge")

    def line(self):
        return "\t".join(self.cols)


class Annotation:
    def __init__(self, path):
        self.features, self.headers, self.fasta = [], [], []
        self.group_directives = 0
        in_fasta = False
        with opener(path) as handle:
            for number, line in enumerate(handle, 1):
                line = line.rstrip("\r\n")
                if in_fasta:
                    self.fasta.append(line)
                elif line == "##FASTA":
                    in_fasta = True
                    self.fasta.append(line)
                elif line == "###":
                    # Adding features changes grouping; omit and report these
                    # positional assertions that all references are resolved.
                    self.group_directives += 1
                elif line.startswith("#"):
                    if line.startswith("##gff-version") and line != "##gff-version 3":
                        raise GFFError(f"line {number}: only GFF3 is supported")
                    self.headers.append(line)
                elif line.strip():
                    try:
                        self.features.append(Feature(line.split("\t"), number))
                    except (ValueError, UnicodeError) as exc:
                        raise GFFError(f"line {number}: {exc}") from exc
        if not self.features:
            raise GFFError("no annotation records found")

    def index(self):
        by_id, children = defaultdict(list), defaultdict(list)
        for f in self.features:
            if f.id:
                by_id[f.id].append(f)
            for p in f.parents:
                children[p].append(f)
        return by_id, children

    def write(self, path):
        with open(path, "w", encoding="utf-8", newline="\n") as handle:
            handle.write("##gff-version 3\n")
            for line in self.headers:
                if not line.startswith("##gff-version"):
                    handle.write(line + "\n")
            for f in self.features:
                handle.write(f.line() + "\n")
            for line in self.fasta:
                handle.write(line + "\n")


def load_roles(extension=None):
    data = json.loads(ROLE_FILE.read_text(encoding="utf-8"))
    roles = {term.lower(): role for role, terms in data["roles"].items() for term in terms}
    if extension:
        extra = json.loads(Path(extension).read_text(encoding="utf-8"))
        if not isinstance(extra, dict):
            raise GFFError("--type-map must be a JSON object: feature_type -> role")
        for term, role in extra.items():
            if not isinstance(role, str) or role not in data["roles"]:
                raise GFFError(f"unsupported role for {term}: {role}")
            if term.lower() in roles and roles[term.lower()] != role:
                raise GFFError("cannot override a bundled role: " + term)
            roles[term.lower()] = role
    return roles


class Audit:
    def __init__(self):
        self.counts, self.examples = Counter(), defaultdict(list)

    def warn(self, code, detail):
        self.counts[code] += 1
        if len(self.examples[code]) < 10:
            self.examples[code].append(detail)

    def report(self):
        return {key: {"count": count, "examples": self.examples[key]}
                for key, count in sorted(self.counts.items())}


def validate(annotation, roles, audit=None, allow_missing=False):
    audit = audit if audit is not None else Audit()
    by_id, children = annotation.index()
    for fid, rows in by_id.items():
        signatures = {(f.type, f.seqid, f.strand, frozenset(f.parents)) for f in rows}
        if len(signatures) != 1:
            raise GFFError(f"incompatible records share ID {fid!r} (types/seqid/strand/Parent differ)")
        if len(rows) > 1:
            audit.warn("segmented-id", fid)
    missing = set(children) - set(by_id)
    if missing and not allow_missing:
        raise GFFError("missing Parent IDs: " + ", ".join(sorted(missing)[:10]) +
                       "; repair requires explicit --repair-missing-parents and consistent transcript evidence")
    # Kahn's algorithm avoids recursion limits on long feature graphs.
    edges, degree = defaultdict(set), {fid: 0 for fid in by_id}
    for fid, rows in by_id.items():
        for parent in rows[0].parents:
            if parent in by_id:
                edges[parent].add(fid)
                degree[fid] += 1
    queue = deque(fid for fid, d in degree.items() if not d)
    visited = 0
    while queue:
        parent = queue.popleft()
        visited += 1
        for child in edges[parent]:
            degree[child] -= 1
            if degree[child] == 0:
                queue.append(child)
    if visited != len(by_id):
        raise GFFError("Parent cycle detected: " + ", ".join(k for k, v in degree.items() if v)[:300])
    for f in annotation.features:
        role = roles.get(f.type.lower())
        if role == "cds" and f.phase == ".":
            raise GFFError(f"CDS requires phase 0, 1 or 2 (line {f.line_number})")
        if role is None:
            audit.warn("unknown-feature-type", f.type)
        if role == "transcript" and not f.parents:
            audit.warn("transcript-without-parent", f.id or f"line {f.line_number}")
        for parent in f.parents:
            if parent in by_id:
                p = by_id[parent][0]
                if roles.get(p.type.lower()) is None:
                    audit.warn("unknown-parent-type", f"{f.id or f.type} -> {parent} ({p.type})")
                if p.seqid != f.seqid:
                    audit.warn("cross-sequence-parent", f"{f.id or f.type} -> {parent}")
                elif not any(p.start <= f.start <= f.end <= p.end for p in by_id[parent]):
                    audit.warn("child-outside-parent", f"{f.id or f.type} -> {parent}")
    return by_id, children


def merge(blocks):
    result = []
    for start, end in sorted(blocks):
        if result and start <= result[-1][1] + 1:
            result[-1][1] = max(result[-1][1], end)
        else:
            result.append([start, end])
    return [tuple(block) for block in result]


def subtract(blocks, holes):
    result = []
    for start, end in blocks:
        pieces = [(start, end)]
        for hs, he in holes:
            new = []
            for s, e in pieces:
                if he < s or hs > e:
                    new.append((s, e))
                else:
                    if s < hs:
                        new.append((s, hs - 1))
                    if e > he:
                        new.append((he + 1, e))
            pieces = new
        result.extend(pieces)
    return result


def normalize(annotation, roles, audit, changes, repair=False):
    by_id, children = validate(annotation, roles, allow_missing=repair)
    used = set(by_id)
    stats = Counter()

    def add(template, kind, start, end, wanted, parent, evidence):
        fid, suffix = wanted, 1
        while fid in used:
            fid = f"{wanted}.cf{suffix}"
            suffix += 1
        used.add(fid)
        attrs = "ID=" + encode(fid)
        if parent:
            attrs += ";Parent=" + encode(parent)
        attrs += ";chromaflow_evidence=" + encode(evidence)
        f = Feature([template.seqid, "ChromaFlow", kind, str(start), str(end), ".",
                     template.strand, ".", attrs])
        annotation.features.append(f)
        stats[kind + "_added"] += 1
        changes.write("\t".join(["add", encode(fid), kind, encode(parent), str(start), str(end), evidence]) + "\n")

    # Recover only explicitly referenced missing genes. Original edges remain
    # untouched, including multiple transcripts sharing the same parent.
    for parent in sorted(set(children) - set(by_id)):
        kids = children[parent]
        if (not repair or any(roles.get(f.type.lower()) != "transcript" for f in kids)
                or len({(f.seqid, f.strand) for f in kids}) != 1
                or any(len(f.parents) != 1 or not f.id for f in kids)):
            raise GFFError(f"cannot safely recover missing Parent {parent!r}")
        add(kids[0], "gene", min(f.start for f in kids), max(f.end for f in kids),
            parent, "", "shared-explicit-Parent;transcript-span")
    for tid, segments in by_id.items():
        tx = segments[0]
        if roles.get(tx.type.lower()) != "transcript":
            continue
        if len(segments) != 1:
            audit.warn("inference-skipped-segmented-transcript", tid)
            continue
        kids = children.get(tid, [])
        blocks = [f for f in kids if roles.get(f.type.lower()) in {"exon", "cds", "utr", "utr5", "utr3"}]
        if not blocks:
            audit.warn("transcript-without-block-evidence", tid)
            continue
        if any(f.seqid != tx.seqid or f.start < tx.start or f.end > tx.end or
               (f.strand in {"+", "-"} and tx.strand in {"+", "-"} and f.strand != tx.strand)
               for f in blocks):
            audit.warn("inference-skipped-special-model", tid)
            continue
        exons = [f for f in blocks if roles.get(f.type.lower()) == "exon"]
        cds = [f for f in blocks if roles.get(f.type.lower()) == "cds"]
        utrs = [f for f in blocks if roles.get(f.type.lower()) in {"utr", "utr5", "utr3"}]
        if not exons:
            for i, (start, end) in enumerate(merge((f.start, f.end) for f in cds + utrs), 1):
                add(tx, "exon", start, end, f"{tid}.exon{i}", tid, "CDS-or-UTR-union")
            continue
        if not cds:
            continue
        if tx.strand not in {"+", "-"} or any(f.strand != tx.strand for f in exons + cds):
            audit.warn("utr-inference-skipped-unknown-strand", tid)
            continue
        exon_blocks = merge((f.start, f.end) for f in exons)
        cds_blocks = merge((f.start, f.end) for f in cds)
        if subtract(cds_blocks, exon_blocks):
            audit.warn("utr-inference-skipped-CDS-outside-exons", tid)
            continue
        remaining = subtract(exon_blocks, merge((f.start, f.end) for f in cds + utrs))
        for i, (start, end) in enumerate(remaining, 1):
            if end < cds_blocks[0][0]:
                kind = "five_prime_UTR" if tx.strand == "+" else "three_prime_UTR"
            elif start > cds_blocks[-1][1]:
                kind = "three_prime_UTR" if tx.strand == "+" else "five_prime_UTR"
            else:
                audit.warn("internal-noncoding-gap-unclassified", f"{tid}:{start}-{end}")
                continue
            add(tx, kind, start, end, f"{tid}.{kind}{i}", tid, "existing-exon-minus-CDS;terminal;known-strand")
    return dict(stats)


def check_preserved(original, result):
    # Counters preserve repeated/anonymous records as well as named features.
    remaining = Counter(f.cols for f in result.features)
    for f in original.features:
        if remaining[f.cols] == 0:
            raise GFFError(f"engine removed/changed source record at line {f.line_number}: {f.id or f.type}; "
                           "use --engine builtin to preserve identity, Parent, coordinates and attributes")
        remaining[f.cols] -= 1
    if original.fasta != result.fasta:
        raise GFFError("engine changed or removed embedded FASTA")


def tool_version(command):
    proc = subprocess.run(command, capture_output=True, text=True, timeout=30)
    if proc.returncode:
        raise GFFError("could not determine tool version: " + " ".join(command))
    return (proc.stdout + proc.stderr).strip()[:2000]


def same_file(a, b):
    a, b = Path(a), Path(b)
    return a.resolve() == b.resolve() or (b.exists() and os.path.samefile(a, b))


def cached(output, metadata, key):
    try:
        data = json.loads(Path(metadata).read_text(encoding="utf-8"))
        return data.get("key") == key and data.get("output_sha256") == digest(output)
    except (OSError, ValueError):
        return False


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def run(args):
    src, dst = Path(args.input).resolve(), Path(args.output).absolute()
    if any(same_file(src, p) for p in (dst, Path(str(dst) + ".report.json"), Path(str(dst) + ".changes.tsv"))):
        raise GFFError("input and output must be different files (including symlinks/hardlinks)")
    roles = load_roles(args.type_map)
    engine = "builtin" if args.engine == "auto" else args.engine
    version = VERSION if engine == "builtin" else tool_version(
        ["gffread", "--version"] if engine == "gffread" else ["agat_convert_sp_gxf2gxf.pl", "--version"])
    if args.repair_missing_parents and engine != "builtin":
        raise GFFError("--repair-missing-parents is supported only by builtin")
    key = {"input_sha256": digest(src), "engine": engine, "engine_version": version,
           "policy": POLICY, "implementation_sha256": digest(__file__),
           "role_table_version": json.loads(ROLE_FILE.read_text(encoding="utf-8"))["version"],
           "roles_sha256": digest(ROLE_FILE), "repair_missing_parents": args.repair_missing_parents,
           "type_map_sha256": digest(args.type_map) if args.type_map else None}
    report_path = Path(str(dst) + ".report.json")
    changes_path = Path(str(dst) + ".changes.tsv")
    if not args.force and cached(dst, report_path, key):
        data = json.loads(report_path.read_text(encoding="utf-8"))
        if changes_path.is_file() and data.get("changes_sha256") == digest(changes_path):
            print("[SKIP] validated cache (input, rules, engine, output and changes hashes match)")
            print("[REPORT]", report_path)
            return
    original = Annotation(src)
    audit = Audit()
    validate(original, roles, allow_missing=args.repair_missing_parents)
    dst.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".chromaflow-gff-", dir=dst.parent) as folder:
        tmp = Path(folder)
        candidate = tmp / "normalized.gff3"
        stats = {}
        with (tmp / "changes.tsv").open("w", encoding="utf-8", newline="\n") as changes:
            changes.write("action\tid\ttype\tparent\tstart\tend\tevidence\n")
            if engine == "builtin":
                model = copy.copy(original)
                model.features = list(original.features)
                stats = normalize(model, roles, audit, changes, args.repair_missing_parents)
                model.write(candidate)
                del model
            else:
                command = (["gffread", str(src), "--keep-genes", "-F", "-O", "-o", str(candidate)]
                           if engine == "gffread" else
                           ["agat_convert_sp_gxf2gxf.pl", "-g", str(src), "-o", str(candidate)])
                subprocess.run(command, check=True)
            result = Annotation(candidate)
            validate(result, roles, audit)
            check_preserved(original, result)
            if engine != "builtin":
                old = Counter(f.cols for f in original.features)
                for f in result.features:
                    if old[f.cols]:
                        old[f.cols] -= 1
                    else:
                        stats[f.type + "_added"] = stats.get(f.type + "_added", 0) + 1
                        changes.write("\t".join(["add", encode(f.id), f.type, ",".join(map(encode, f.parents)),
                                                  str(f.start), str(f.end), "external-engine"]) + "\n")
        if digest(src) != key["input_sha256"]:
            raise GFFError("input changed during normalization; output was not published")
        report = {"status": "validated", "scope": "GFF3 syntax, Parent graph and source-record preservation; not downstream validation",
                  "chromaflow_version": VERSION, "key": key, "input": str(src), "output": str(dst),
                  "input_records": len(original.features), "output_records": len(result.features),
                  "original_records_preserved": True, "changes": stats,
                  "warnings": audit.report(), "removed_group_directives": original.group_directives,
                  "feature_counts": dict(Counter(f.type for f in result.features)),
                  "output_sha256": digest(candidate), "changes_sha256": digest(tmp / "changes.tsv")}
        write_json(tmp / "report.json", report)
        os.replace(tmp / "changes.tsv", changes_path)
        os.replace(tmp / "report.json", report_path)
        os.replace(candidate, dst)
    print("[OK] original IDs, Parent edges, coordinates, types, phases and attributes preserved")
    print("[ADDED]", json.dumps(stats, ensure_ascii=False))
    for code, entry in report["warnings"].items():
        print(f"[WARN] {code}: {entry['count']}; examples: {entry['examples']}")
    print("[REPORT]", report_path)
    print("[CHANGES]", changes_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("--engine", choices=["auto", "builtin", "agat", "gffread"], default="builtin")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--repair-missing-parents", action="store_true")
    parser.add_argument("--type-map", help="local JSON mapping: feature_type -> role")
    args = parser.parse_args()
    try:
        run(args)
    except (ValueError, OSError, UnicodeError, EOFError, zipfile.BadZipFile, subprocess.SubprocessError) as exc:
        print("[ERROR] " + str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
