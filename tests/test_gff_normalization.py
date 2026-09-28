"""Structural regressions, publication safety and downstream identity checks."""
import gzip
import bz2
import io
import json
import lzma
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import zipfile

STEPS = Path(__file__).resolve().parents[1] / "lib/steps"
sys.path.insert(0, str(STEPS))
from normalize_gff import Annotation, Audit, GFFError, check_preserved, load_roles, normalize, validate
from prepare_gtf import check_gtf, check_source, expected_assignments, expected_exons


def row(kind, fid="", parent="", start=1, end=100, strand="+", seqid="chr1", extra=""):
    attrs = (["ID=" + fid] if fid else []) + (["Parent=" + parent] if parent else [])
    if extra:
        attrs.append(extra)
    return "\t".join([seqid, "fixture", kind, str(start), str(end), ".", strand,
                      "0" if kind in {"CDS", "SO:0000316"} else ".", ";".join(attrs) or "."])


def basic(gene="gene", strand="+"):
    return [row(gene, "g", strand=strand), row("mRNA", "t", "g", strand=strand),
            row("CDS", "c", "t", 20, 80, strand=strand)]


class GFFTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="cf-gff-test-")
        self.root = Path(self.tmp.name)
        self.src, self.dst = self.root / "input.gff3", self.root / "out.gff3"
        self.roles = load_roles()

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, lines):
        self.src.write_text("##gff-version 3\n" + "\n".join(lines) + "\n", encoding="utf-8")

    def model(self, lines, repair=False):
        self.write(lines)
        old, new = Annotation(self.src), Annotation(self.src)
        audit, changes = Audit(), io.StringIO()
        normalize(new, self.roles, audit, changes, repair)
        validate(new, self.roles, audit)
        check_preserved(old, new)
        return new, audit

    def cli(self, *extra, success=True):
        p = subprocess.run([sys.executable, str(STEPS / "normalize_gff.py"), str(self.src), str(self.dst), *extra],
                           capture_output=True, text=True, encoding="utf-8")
        self.assertEqual(p.returncode == 0, success, p.stdout + p.stderr)
        return p

    def test_te_gene_relationships_all_preserved(self):
        # 810 independent synthetic TE models, not the user's unavailable dataset.
        lines = []
        for i in range(810):
            lines += [row("transposable_element_gene", f"TE{i}"), row("mRNA", f"t{i}", f"TE{i}"),
                      row("CDS", parent=f"t{i}", start=20, end=80)]
        model, _ = self.model(lines)
        self.assertEqual(sum(f.type == "gene" for f in model.features), 0)
        self.assertEqual(sum(f.type == "exon" for f in model.features), 810)
        self.assertEqual(len(expected_assignments(model, self.roles)), 810)

    def test_unknown_existing_parent_never_replaced(self):
        model, audit = self.model(basic("provider_gene_model"))
        self.assertEqual(model.features[1].parents, ("g",))
        self.assertEqual(sum(f.type == "gene" for f in model.features), 0)
        self.assertIn("unknown-parent-type", audit.counts)
        with self.assertRaisesRegex(GFFError, "cannot interpret"):
            expected_assignments(model, self.roles)

    def test_missing_parent_is_error_by_default(self):
        with self.assertRaisesRegex(GFFError, "missing Parent"):
            self.model(basic()[1:])

    def test_explicit_repair_keeps_shared_missing_id(self):
        model, _ = self.model([row("mRNA", "a", "g"), row("mRNA", "b", "g", end=200)], repair=True)
        genes = [f for f in model.features if f.type == "gene"]
        self.assertEqual([(f.id, f.start, f.end) for f in genes], [("g", 1, 200)])
        self.assertEqual([f.parents for f in model.features[:2]], [("g",), ("g",)])

    def test_repair_refuses_inconsistent_evidence(self):
        for bad in (row("mRNA", "b", "g", seqid="other"), row("exon", "b", "g"),
                    row("mRNA", "b", "g", strand="-")):
            with self.subTest(bad=bad), self.assertRaisesRegex(GFFError, "cannot safely recover"):
                self.model([row("mRNA", "a", "g"), bad], repair=True)

    def test_original_spans_and_source_attributes_preserved(self):
        model, _ = self.model([row("gene", "g", end=1000), row("mRNA", "t", "g", end=900),
                               row("CDS", "c", "t", 100, 600, extra="Name=A%3BB;Dbxref=X:1;Note=some%20text")])
        self.assertEqual(model.features[0].end, 1000)
        self.assertEqual(model.features[1].end, 900)
        self.assertEqual(model.features[2].phase, "0")

    def test_no_fabricated_exon_without_evidence(self):
        model, audit = self.model(basic()[:2])
        self.assertEqual(len(model.features), 2)
        self.assertIn("transcript-without-block-evidence", audit.counts)

    def test_unknown_strand_never_infers_directional_utr(self):
        for strand in (".", "?"):
            model, _ = self.model(basic(strand=strand) + [row("exon", "e", "t", strand=strand)])
            self.assertFalse(any("UTR" in f.type for f in model.features))

    def test_terminal_utr_orientation_and_internal_gap(self):
        for strand in ("+", "-"):
            lines = basic(strand=strand)[:2] + [row("exon", "e", "t", strand=strand),
                     row("CDS", "c", "t", 20, 30, strand=strand), row("CDS", "c", "t", 70, 80, strand=strand)]
            model, audit = self.model(lines)
            utr = [(f.type, f.start, f.end) for f in model.features if "UTR" in f.type]
            self.assertEqual(len(utr), 2)
            self.assertIn(("five_prime_UTR", 1, 19) if strand == "+" else ("three_prime_UTR", 1, 19), utr)
            self.assertIn("internal-noncoding-gap-unclassified", audit.counts)

    def test_partial_utr_only_adds_uncovered_region(self):
        model, _ = self.model(basic() + [row("exon", "e", "t"), row("five_prime_UTR", "u", "t", 1, 10)])
        self.assertEqual([(f.start, f.end) for f in model.features if f.source == "ChromaFlow" and f.type == "five_prime_UTR"], [(11, 19)])

    def test_generic_utr_stays_untyped_and_not_duplicated(self):
        model, _ = self.model(basic() + [row("exon", "e", "t"), row("UTR", "u", "t", 1, 19)])
        self.assertEqual(sum(f.type == "UTR" for f in model.features), 1)
        self.assertEqual(sum(f.type == "five_prime_UTR" for f in model.features), 0)

    def test_ncRNA_needs_no_CDS_or_UTR(self):
        model, _ = self.model([row("ncRNA_gene", "g"), row("lnc_RNA", "t", "g"), row("exon", "e", "t")])
        self.assertEqual(len(model.features), 3)

    def test_SO_accessions_and_custom_role_mapping(self):
        model, _ = self.model([row("SO:0000704", "g"), row("SO:0000234", "t", "g"), row("SO:0000316", "c", "t")])
        self.assertEqual(model.features[-1].type, "exon")
        mapping = self.root / "types.json"
        mapping.write_text('{"custom_gene":"gene"}')
        self.assertEqual(load_roles(mapping)["custom_gene"], "gene")
        mapping.write_text('{"gene":"transcript"}')
        with self.assertRaisesRegex(GFFError, "cannot override"):
            load_roles(mapping)

    def test_multi_parent_exon_and_direct_gene_CDS_preserved(self):
        direct = [row("gene", "g"), row("gene", "h"), row("CDS", "c", "g,h")]
        model, _ = self.model(direct)
        self.assertEqual(len(model.features), 3)
        with self.assertRaisesRegex(GFFError, "direct gene-to"):
            expected_assignments(model, self.roles)
        model, _ = self.model(direct + [row("mRNA", "t", "g"), row("mRNA", "u", "h"), row("exon", "e", "t,u")])
        self.assertEqual(model.features[2].parents, ("g", "h"))
        self.assertEqual(model.features[-1].parents, ("t", "u"))

    def test_legal_segmented_CDS_id(self):
        model, audit = self.model(basic()[:2] + [row("CDS", "c", "t", 10, 30), row("CDS", "c", "t", 60, 90)])
        self.assertEqual(sum(f.id == "c" for f in model.features), 2)
        self.assertIn("segmented-id", audit.counts)

    def test_incompatible_duplicate_id_fails(self):
        with self.assertRaisesRegex(GFFError, "incompatible records"):
            self.model(basic() + [row("exon", "g", "t")])

    def test_inferred_id_never_collides(self):
        model, _ = self.model(basic() + [row("repeat_region", "t.exon1")])
        self.assertIn("t.exon1.cf1", [f.id for f in model.features])

    def test_malformed_input_is_fatal(self):
        cases = [row("CDS", "c", "t", start="bad"), row("CDS", "c", "t", start=0),
                 row("CDS", "c", "t", start=200), row("CDS", "c", "t", strand="x"),
                 row("CDS", "c", "t") + "\textra", row("CDS", "c", "t", extra="ID=again"),
                 row("CDS", "c", "t", extra="Note=bad%q"), "chr1\tbad\tgene"]
        for case in cases:
            with self.subTest(case=case), self.assertRaisesRegex(GFFError, "line 4"):
                self.model(basic()[:2] + [case])

    def test_cycles_fatal(self):
        with self.assertRaisesRegex(GFFError, "cycle"):
            self.model([row("gene", "g", "t"), row("mRNA", "t", "g")])

    def test_escaped_commas_are_not_parent_separators(self):
        model, _ = self.model([row("gene", "g%2C1"), row("mRNA", "t%3B1", "g%2C1"), row("CDS", "c", "t%3B1")])
        self.assertEqual(model.features[1].parents, ("g,1",))
        self.assertEqual(model.features[-1].parents, ("t;1",))

    def test_embedded_fasta_and_comments_preserved(self):
        self.write(["# provider information", *basic(), "###", "##FASTA", ">chr1", "ACGT"])
        self.cli()
        text = self.dst.read_text()
        self.assertIn("# provider information", text)
        self.assertTrue(text.endswith("##FASTA\n>chr1\nACGT\n"))
        self.assertNotIn("###", text)

    def test_gzip_magic_without_extension(self):
        self.write(basic())
        self.src.write_bytes(gzip.compress(self.src.read_bytes()))
        self.cli()
        self.assertEqual(len(Annotation(self.dst).features), 4)

    def test_normalization_is_idempotent(self):
        for lines in (basic(), basic() + [row("exon", "e", "t")]):
            self.write(lines)
            self.cli("--force")
            first = self.dst.read_bytes()
            self.src.write_bytes(first)
            self.cli("--force")
            self.assertEqual(first, self.dst.read_bytes())

    def test_cache_uses_contents_and_rejects_tampering(self):
        self.write(basic())
        self.cli()
        self.assertIn("[SKIP]", self.cli().stdout)
        self.dst.write_text("corrupted")
        self.assertNotIn("[SKIP]", self.cli().stdout)
        self.write(basic("transposable_element_gene"))
        self.assertNotIn("[SKIP]", self.cli().stdout)

    def test_failed_rebuild_preserves_previous_product_and_report(self):
        self.write(basic())
        self.cli()
        previous = {p: p.read_bytes() for p in self.root.glob("out.gff3*")}
        self.write(basic()[1:])
        self.cli("--force", success=False)
        for p, data in previous.items():
            self.assertEqual(data, p.read_bytes())

    def test_same_file_and_hardlink_rejected(self):
        self.write(basic())
        original = self.src.read_bytes()
        os.link(self.src, self.dst)
        self.cli("--force", success=False)
        self.assertEqual(original, self.src.read_bytes())

    def test_all_engines_must_preserve_original_record(self):
        self.write(basic())
        original = Annotation(self.src)
        self.write([row("gene", "g"), row("mRNA", "t", "g", start=2), basic()[2]])
        with self.assertRaisesRegex(GFFError, "removed/changed"):
            check_preserved(original, Annotation(self.src))

    def test_gtf_identity_checker_detects_reassigned_or_missing_transcript(self):
        self.dst.write_text('chr1\tgffread\texon\t1\t100\t.\t+\t.\tgene_id "g"; transcript_id "t";\n')
        self.assertTrue(check_gtf(self.dst, {"t": {"g"}})["transcript_gene_assignments_preserved"])
        with self.assertRaisesRegex(GFFError, "changed gene"):
            check_gtf(self.dst, {"t": {"TE"}})
        with self.assertRaisesRegex(GFFError, "lost transcript"):
            check_gtf(self.dst, {"t": {"g"}, "missing": {"g"}})

    def test_gtf_without_final_semicolon(self):
        self.dst.write_text('chr1\tgffread\texon\t1\t100\t.\t+\t.\ttranscript_id "t"; gene_id "g"\n')
        self.assertEqual(check_gtf(self.dst, {"t": {"g"}})["transcripts"], 1)

    def test_other_compressions_are_read_by_magic(self):
        for compress in (bz2.compress, lzma.compress):
            self.write(basic())
            self.src.write_bytes(compress(self.src.read_bytes()))
            self.cli("--force")
        self.write(basic())
        raw = self.src.read_bytes()
        with zipfile.ZipFile(self.src, "w") as archive:
            archive.writestr("genes.gff3", raw)
        self.cli("--force")

    def test_gtf_rejects_unsupported_exon_invention_or_span_change(self):
        model, _ = self.model(basic()[:2])
        with self.assertRaisesRegex(GFFError, "no exon/CDS/UTR evidence"):
            expected_exons(model, self.roles, {"t": {"g"}})
        self.dst.write_text('chr1\tgffread\texon\t1\t100\t.\t+\t.\ttranscript_id "t"; gene_id "g"\n')
        with self.assertRaisesRegex(GFFError, "exon spans differ"):
            check_gtf(self.dst, {"t": {"g"}}, {"t": ("chr1", "+", [(20, 80)])})

    def test_source_and_rules_changes_invalidate_normalization(self):
        self.write(basic())
        self.cli()
        check_source(self.dst, self.src)
        self.write(basic("transposable_element_gene"))
        with self.assertRaisesRegex(GFFError, "does not match"):
            check_source(self.dst, self.src)
        report_path = Path(str(self.dst) + ".report.json")
        report = json.loads(report_path.read_text())
        report["key"]["implementation_sha256"] = "old rules"
        report_path.write_text(json.dumps(report))
        with self.assertRaisesRegex(GFFError, "rules/type map changed"):
            check_source(self.dst, None)

    def test_output_sidecar_cannot_overwrite_input(self):
        self.src = Path(str(self.dst) + ".report.json")
        self.write(basic())
        previous = self.src.read_bytes()
        self.cli(success=False)
        self.assertEqual(self.src.read_bytes(), previous)


if __name__ == "__main__":
    unittest.main()
