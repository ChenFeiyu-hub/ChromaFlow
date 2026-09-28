#!/usr/bin/env python3
"""Content-addressed GFF -> GTF conversion with transcript/gene identity checks."""
import argparse
from collections import defaultdict
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

from normalize_gff import (Annotation, Audit, GFFError, ROLE_FILE, VERSION, cached,
                           digest, load_roles, merge, same_file, tool_version, validate, write_json)


def expected_assignments(model, roles):
    by_id, children = validate(model, roles)
    expected = {}
    for f in model.features:
        role = roles.get(f.type.lower())
        if role == "transcript":
            if not f.id or not f.parents:
                raise GFFError(f"GTF export requires explicit transcript ID and gene Parent (line {f.line_number})")
            if any(roles.get(by_id[p][0].type.lower()) != "gene" for p in f.parents):
                raise GFFError(f"cannot interpret gene parents for transcript {f.id!r}; "
                               "GFF was preserved; supply an evidence-based GFF_TYPE_MAP or a separate compatible annotation")
            expected[f.id] = set(f.parents)
        elif role == "gene" and any(roles.get(k.type.lower()) in {"cds", "exon"} for k in children.get(f.id, [])):
            # Legal GFF3, but converter-specific transcript synthesis is not a
            # universally reversible export. Do not silently make up a mapping.
            raise GFFError(f"direct gene-to-CDS/exon model {f.id!r} is preserved in GFF; "
                           "automatic GTF export of this structure is not supported by the identity checker")
    # An unknown parent with exon/CDS children is a potential unrecognized
    # transcript. Dropping it silently would defeat the preservation check.
    for parent, kids in children.items():
        if roles.get(by_id[parent][0].type.lower()) is None and any(
                roles.get(k.type.lower()) in {"cds", "exon"} for k in kids):
            raise GFFError(f"unknown model type at {parent!r}; provide GFF_TYPE_MAP before GTF export")
    if not expected:
        raise GFFError("no supported transcript-to-gene assignments for GTF export")
    return expected


def expected_exons(model, roles, expected):
    _, children = model.index()
    blocks = {}
    for tid in expected:
        kids = children.get(tid, [])
        exons = [f for f in kids if roles.get(f.type.lower()) == "exon"]
        evidence = exons or [f for f in kids if roles.get(f.type.lower()) in {"cds", "utr", "utr5", "utr3"}]
        if not evidence:
            raise GFFError(f"GTF export has no exon/CDS/UTR evidence for {tid!r}; "
                           "refusing a converter-invented whole-transcript exon")
        if len({(f.seqid, f.strand) for f in evidence}) != 1:
            raise GFFError(f"GTF export of a cross-sequence/strand model is not supported: {tid!r}")
        blocks[tid] = (evidence[0].seqid, evidence[0].strand, merge((f.start, f.end) for f in evidence))
    return blocks


def check_gtf(path, expected, expected_blocks=None):
    observed = defaultdict(set)
    observed_blocks = defaultdict(list)
    records = 0
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) != 9:
                raise GFFError(f"invalid GTF row at line {number}")
            attrs = dict(re.findall(r'(\w+)\s+"([^"\n]*)"\s*(?:;|$)', cols[8]))
            tid, gid = attrs.get("transcript_id"), attrs.get("gene_id")
            if not tid or not gid:
                raise GFFError(f"GTF line {number} lacks transcript_id/gene_id")
            if tid not in expected or gid not in expected[tid]:
                raise GFFError(f"GTF changed gene assignment: {tid!r} -> {gid!r}")
            observed[tid].add(gid)
            if expected_blocks is not None and cols[2] == "exon":
                seqid, strand, _ = expected_blocks[tid]
                if cols[0] != seqid or cols[6] != strand:
                    raise GFFError(f"GTF changed exon sequence/strand for {tid!r}")
                start, end = int(cols[3]), int(cols[4])
                if start < 1 or end < start:
                    raise GFFError(f"invalid GTF exon coordinates at line {number}")
                observed_blocks[tid].append((start, end))
            records += 1
    missing = [(tid, sorted(genes - observed[tid])) for tid, genes in expected.items() if genes != observed[tid]]
    if missing:
        raise GFFError(f"GTF lost transcript/gene assignments: {missing[:10]} (total {len(missing)})")
    if expected_blocks is not None:
        bad = [tid for tid, (_, _, blocks) in expected_blocks.items() if merge(observed_blocks[tid]) != blocks]
        if bad:
            raise GFFError(f"GTF exon spans differ from source evidence: {bad[:10]} (total {len(bad)})")
    return {"records": records, "transcripts": len(observed),
            "genes": len(set().union(*observed.values())), "transcript_gene_assignments_preserved": True,
            "exon_spans_preserved": True if expected_blocks is not None else None}


def check_source(gff, source, type_map=None):
    report = Path(str(gff) + ".report.json")
    if not report.exists():
        return
    data = json.loads(report.read_text(encoding="utf-8"))
    key = data.get("key", {})
    if (key.get("implementation_sha256") != digest(Path(__file__).with_name("normalize_gff.py"))
            or key.get("roles_sha256") != digest(ROLE_FILE)
            or key.get("type_map_sha256") != (digest(type_map) if type_map else None)):
        raise GFFError("normalization rules/type map changed; rerun normalize-gff")
    if data.get("output_sha256") != digest(gff):
        raise GFFError("normalized GFF differs from its report; rerun normalize-gff")
    if source and not same_file(gff, source) and data.get("key", {}).get("input_sha256") != digest(source):
        raise GFFError("normalized GFF does not match current GFF3 source; rerun normalize-gff")


def run(args):
    src, dst = Path(args.input).resolve(), Path(args.output).absolute()
    check_source(src, args.source, args.type_map)
    if args.check_source_only:
        return
    if same_file(src, dst) or same_file(src, Path(str(dst) + ".report.json")):
        raise GFFError("GFF input and GTF output must be different files")
    roles = load_roles(args.type_map)
    key = {"input_sha256": digest(src), "gffread_version": tool_version(["gffread", "--version"]),
           "implementation_sha256": digest(__file__), "parser_sha256": digest(Path(__file__).with_name("normalize_gff.py")),
           "roles_sha256": digest(ROLE_FILE), "type_map_sha256": digest(args.type_map) if args.type_map else None}
    report_path = Path(str(dst) + ".report.json")
    if not args.force and cached(dst, report_path, key):
        print("[SKIP] GTF validated cache (GFF, converter, rules and output hashes match)")
        return
    annotation = Annotation(src)
    audit = Audit()
    validate(annotation, roles, audit)
    expected = expected_assignments(annotation, roles)
    blocks = expected_exons(annotation, roles, expected)
    dst.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".chromaflow-gtf-", dir=dst.parent) as folder:
        folder = Path(folder)
        candidate = folder / "genes.gtf"
        subprocess.run(["gffread", str(src), "-T", "-o", str(candidate)], check=True)
        stats = check_gtf(candidate, expected, blocks)
        if key["input_sha256"] != digest(src):
            raise GFFError("GFF changed during conversion; output was not published")
        report = {"status": "validated", "scope": "GTF transcript/gene identity and exon spans; not TxDb/TSS or sequence validation",
                  "chromaflow_version": VERSION, "key": key, "input": str(src), "output": str(dst),
                  "output_sha256": digest(candidate), "statistics": stats, "warnings": audit.report()}
        write_json(folder / "report.json", report)
        os.replace(folder / "report.json", report_path)
        os.replace(candidate, dst)
    print("[OK] GTF transcript/gene identity preserved:", json.dumps(stats))
    print("[REPORT]", report_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--type-map")
    parser.add_argument("--source", help="original configured GFF3 source, to detect stale normalization")
    parser.add_argument("--check-source-only", action="store_true")
    args = parser.parse_args()
    try:
        run(args)
    except (ValueError, OSError, UnicodeError, EOFError, subprocess.SubprocessError) as exc:
        print("[ERROR] " + str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
