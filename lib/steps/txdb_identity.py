#!/usr/bin/env python3
"""Validate TxDb identity through explicit GFF attributes, never suffix guesses.

Works on exports from the installed R importer. Keeps original GFF IDs as the
canonical identity while allowing the importer's public names to differ.
"""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import sys

from normalize_gff import Annotation, GFFError, decode, digest, load_roles, write_json
from prepare_gtf import expected_assignments

POLICY = "explicit-gff-attributes-v1"


def read_tsv(path, required):
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not set(required).issubset(reader.fieldnames or []):
            raise GFFError(f"missing columns in {path}: {required}")
        rows = list(reader)
        if any(None in row or any(row[k] is None for k in required) for row in rows):
            raise GFFError(f"malformed TSV record in {path}")
        return rows


def write_tsv(path, fields, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def alias_index(model, roles, role):
    fields = ("ID", "Name", "transcript_id") if role == "transcript" else ("ID", "Name", "gene_id")
    index, evidence = defaultdict(set), defaultdict(set)
    for f in model.features:
        if roles.get(f.type.lower()) != role or not f.id:
            continue
        for key in fields:
            value = decode(f.attrs.get(key, ""))
            if value:
                index[value].add(f.id)
                evidence[(value, f.id)].add(key)
    return index, evidence


def resolve(value, index, kind):
    candidates = index.get(value, set())
    if len(candidates) != 1:
        reason = "ambiguous" if candidates else "unmapped"
        raise GFFError(f"{reason} {kind} alias {value!r}" +
                       (": " + ", ".join(sorted(candidates)[:8]) if candidates else ""))
    return next(iter(candidates))


def check_exports(model, roles, export_dir, output_dir=None, gff_path=None):
    export_dir = Path(export_dir)
    output_dir = Path(output_dir) if output_dir else export_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    # A failed retry must not leave an old successful summary as current.
    write_json(output_dir / "TxDb-identity-summary.json", {"status": "not_validated", "policy": POLICY})
    expected = expected_assignments(model, roles)
    by_id, _ = model.index()
    tx_aliases, tx_evidence = alias_index(model, roles, "transcript")
    gene_aliases, gene_evidence = alias_index(model, roles, "gene")
    observed_path = export_dir / "TxDb-transcripts.tsv"
    observed = read_tsv(observed_path, ("tx_id", "seqid", "strand", "tss", "gene_ids"))
    edges_path = export_dir / "TxDb-transcript-genes.tsv"
    edges = defaultdict(list)
    if edges_path.is_file():
        for edge in read_tsv(edges_path, ("txdb_id", "gene_id")):
            edges[edge["txdb_id"]].append(edge["gene_id"])
        keys = [row.get("txdb_id", "") for row in observed]
        if any(not k for k in keys) or len(keys) != len(set(keys)):
            raise GFFError("TxDb internal IDs must be present and unique for the gene-edge table")
        if set(edges) - set(keys):
            raise GFFError("gene-edge table references missing TxDb internal IDs")
    elif any("," in name for name in gene_aliases):
        raise GFFError("legacy comma-joined export cannot safely represent these gene aliases; rerun the patched R exporter")

    problems, tx_map, gene_map = [], [], {}
    seen = Counter()
    renamed_tx, renamed_genes = set(), set()
    matched_associations = 0

    def fail(tid, label, detail=""):
        problems.append({"tx_id": tid, "problem": label, "detail": detail})

    for row in observed:
        public_tx = row["tx_id"]  # Historical export header; this is TxDb tx_name.
        try:
            tid = resolve(public_tx, tx_aliases, "transcript")
        except GFFError as exc:
            fail(public_tx, "transcript alias unresolved", str(exc))
            continue
        seen[tid] += 1
        if seen[tid] > 1:
            fail(tid, "duplicate TxDb transcript", public_tx)
            continue
        if tid not in expected:
            fail(tid, "unexpected TxDb transcript", public_tx)
            continue
        segments = by_id[tid]
        if len(segments) != 1 or segments[0].strand not in {"+", "-"}:
            fail(tid, "unsupported ambiguous transcript span/strand")
            continue
        tx = segments[0]
        wanted_tss = tx.start if tx.strand == "+" else tx.end
        public_genes = (edges.get(row["txdb_id"], []) if edges_path.is_file()
                        else [g for g in row["gene_ids"].split(",") if g])
        canonical_genes = set()
        valid = True
        for public_gene in public_genes:
            try:
                gid = resolve(public_gene, gene_aliases, "gene")
            except GFFError as exc:
                fail(tid, "gene alias unresolved", str(exc))
                valid = False
                continue
            canonical_genes.add(gid)
            gene_map[(gid, public_gene)] = {
                "gff_gene_id": gid, "txdb_gene_id": public_gene,
                "matched_attributes": ",".join(sorted(gene_evidence[(public_gene, gid)]))}
            if gid != public_gene:
                renamed_genes.add(gid)
        if canonical_genes != expected[tid]:
            fail(tid, "gene assignment changed", json.dumps({"expected": sorted(expected[tid]),
                                                            "observed": sorted(canonical_genes)}, ensure_ascii=False))
            valid = False
        try:
            same = row["seqid"] == tx.seqid and row["strand"] == tx.strand and int(row["tss"]) == wanted_tss
            if "tx_start" in row and "tx_end" in row:
                same = same and int(row["tx_start"]) == tx.start and int(row["tx_end"]) == tx.end
        except (TypeError, ValueError):
            same = False
        if not same:
            fail(tid, "sequence/strand/TSS/span changed", json.dumps(row, ensure_ascii=False))
            valid = False
        if tid != public_tx:
            renamed_tx.add(tid)
        tx_map.append({"gff_tx_id": tid, "txdb_tx_name": public_tx,
                       "matched_attributes": ",".join(sorted(tx_evidence[(public_tx, tid)])),
                       "gff_gene_ids": json.dumps(sorted(expected[tid]), ensure_ascii=False),
                       "txdb_gene_ids": json.dumps(sorted(set(public_genes)), ensure_ascii=False),
                       "seqid": tx.seqid, "strand": tx.strand, "tss": wanted_tss,
                       "status": "validated" if valid else "failed"})
        if valid:
            matched_associations += len(expected[tid])
    for tid in sorted(set(expected) - set(seen)):
        fail(tid, "transcript missing after explicit alias mapping")
    # The importer must not split one canonical gene over multiple public IDs.
    public_by_gene = defaultdict(set)
    for gid, public_gene in gene_map:
        public_by_gene[gid].add(public_gene)
    for gid, names in public_by_gene.items():
        if len(names) > 1:
            fail("", "gene identity split across TxDb names", gid + ": " + ", ".join(sorted(names)))
    summary = {
        "status": "failed" if problems else "validated", "policy": POLICY,
        "scope": "explicit GFF aliases, exact gene relationships, sequence/strand/TSS; spans when exported",
        "expected_transcripts": len(expected), "observed_transcripts": len(observed),
        "expected_associations": sum(map(len, expected.values())), "matched_associations": matched_associations,
        "transcripts_with_different_public_name": len(renamed_tx),
        "genes_with_different_public_name": len(renamed_genes),
        "problem_count": len(problems), "problem_counts": dict(Counter(p["problem"] for p in problems)),
        "gff_sha256": digest(gff_path) if gff_path else None,
        "transcript_export_sha256": digest(observed_path),
        "gene_edge_export_sha256": digest(edges_path) if edges_path.is_file() else None,
    }
    write_tsv(output_dir / "TxDb-identity-map.tsv",
              ("gff_tx_id", "txdb_tx_name", "matched_attributes", "gff_gene_ids", "txdb_gene_ids", "seqid", "strand", "tss", "status"), tx_map)
    write_tsv(output_dir / "TxDb-gene-identity-map.tsv",
              ("gff_gene_id", "txdb_gene_id", "matched_attributes"), [gene_map[k] for k in sorted(gene_map)])
    write_tsv(output_dir / "TxDb-problems.tsv", ("tx_id", "problem", "detail"), problems)
    write_json(output_dir / "TxDb-identity-summary.json", summary)
    if problems:
        raise GFFError(f"TxDb identity acceptance failed: {summary['problem_counts']}; see {output_dir / 'TxDb-problems.tsv'}")
    print(f"[OK] TxDb: {len(expected)} transcripts; {matched_associations} gene associations; sequence/strand/TSS preserved")
    print(f"[MAP] explicit aliases: {len(renamed_tx)} transcript names and {len(renamed_genes)} gene names differ from GFF IDs")
    print("[REPORT]", output_dir / "TxDb-identity-summary.json")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gff")
    parser.add_argument("export_dir")
    parser.add_argument("--output-dir", required=True, help="separate directory for replay results")
    parser.add_argument("--type-map")
    args = parser.parse_args()
    try:
        check_exports(Annotation(args.gff), load_roles(args.type_map), args.export_dir, args.output_dir, args.gff)
        print("[NOTE] Existing TSV exports were checked; this command did not rebuild TxDb.")
    except (ValueError, OSError, UnicodeError) as exc:
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        write_json(Path(args.output_dir) / "replay-status.json", {"status": "failed", "error": str(exc)})
        print("[ERROR] " + str(exc), file=sys.stderr)
        return 1
    write_json(Path(args.output_dir) / "replay-status.json", {"status": "validated", "mode": "existing exports only"})
    return 0


if __name__ == "__main__":
    sys.exit(main())
