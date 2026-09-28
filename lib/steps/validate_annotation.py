#!/usr/bin/env python3
"""Read-only acceptance of normalized GFF, GTF and optional installed TxDb."""
import argparse
import csv
from pathlib import Path
import subprocess
import sys

from normalize_gff import Annotation, Audit, GFFError, check_preserved, digest, load_roles, validate, write_json
from prepare_gtf import check_gtf, expected_assignments, expected_exons
from txdb_identity import check_exports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("original")
    parser.add_argument("normalized")
    parser.add_argument("--gtf")
    parser.add_argument("--txdb", action="store_true", help="run the installed R/TxDb stack and check gene links and TSS")
    parser.add_argument("--type-map")
    parser.add_argument("--report-dir", required=True)
    args = parser.parse_args()
    out = Path(args.report_dir).absolute()
    try:
        roles, audit = load_roles(args.type_map), Audit()
        original, normalized = Annotation(args.original), Annotation(args.normalized)
        validate(original, roles)
        validate(normalized, roles, audit)
        check_preserved(original, normalized)
        report = {"original_sha256": digest(args.original), "normalized_sha256": digest(args.normalized),
                  "original_records_preserved": True, "original_records": len(original.features),
                  "normalized_records": len(normalized.features), "warnings": audit.report(),
                  "gtf": "not requested", "txdb": "not requested"}
        expected = expected_assignments(normalized, roles) if args.gtf or args.txdb else {}
        if args.gtf:
            report["gtf"] = check_gtf(args.gtf, expected, expected_exons(normalized, roles, expected))
            report["gtf_sha256"] = digest(args.gtf)
        out.mkdir(parents=True, exist_ok=True)
        if args.txdb:
            by_id, _ = normalized.index()
            expected_path = out / "expected_transcripts.tsv"
            with expected_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["tx_id", "gene_id", "seqid", "strand", "tss"])
                for tid, genes in sorted(expected.items()):
                    segments = by_id[tid]
                    if len(segments) != 1 or segments[0].strand not in {"+", "-"}:
                        raise GFFError(f"TxDb TSS acceptance requires an unambiguous single-strand transcript: {tid}")
                    tx = segments[0]
                    for gene in sorted(genes):
                        writer.writerow([tid, gene, tx.seqid, tx.strand, tx.start if tx.strand == "+" else tx.end])
            subprocess.run(["Rscript", str(Path(__file__).with_name("check_txdb.R")),
                            str(Path(args.normalized).resolve()), str(expected_path), str(out),
                            "export-for-identity-v2"], check=True)
            report["txdb"] = check_exports(normalized, roles, out, gff_path=args.normalized)
        write_json(out / "acceptance.json", report)
        print("[OK] source records/relationships preserved; requested downstream checks passed")
        print("[REPORT]", out / "acceptance.json")
    except (ValueError, OSError, UnicodeError, subprocess.SubprocessError) as exc:
        # Never leave a stale success report after a failed repeat of acceptance.
        out.mkdir(parents=True, exist_ok=True)
        write_json(out / "acceptance.json", {"status": "failed", "error": str(exc)})
        print("[ERROR] " + str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
