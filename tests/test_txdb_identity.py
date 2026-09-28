"""Exercise explicit name mappings without requiring Bioconductor locally."""
import contextlib
import csv
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib/steps"))
from normalize_gff import Annotation, GFFError, load_roles
from txdb_identity import check_exports, write_tsv
from test_gff_normalization import row


class TxDbIdentityTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="cf-txdb-")
        self.root = Path(self.tmp.name)
        self.gff = self.root / "input.gff3"
        self.out = self.root / "audit"
        self.roles = load_roles()

    def tearDown(self):
        self.tmp.cleanup()

    def model(self, records):
        self.gff.write_text("##gff-version 3\n" + "\n".join(records) + "\n", encoding="utf-8")
        return Annotation(self.gff)

    def observed(self, rows, edges=None):
        # row = public transcript name, sequence, strand, TSS, public gene names
        data = [dict(zip(("tx_id", "seqid", "strand", "tss", "gene_ids"), r)) for r in rows]
        fields = ["tx_id", "seqid", "strand", "tss", "gene_ids"]
        if edges is not None:
            fields += ["txdb_id"]
            for i, record in enumerate(data, 1):
                record["txdb_id"] = str(i)
            write_tsv(self.root / "TxDb-transcript-genes.tsv", ("txdb_id", "gene_id"),
                      [{"txdb_id": str(key), "gene_id": value} for key, value in edges])
        write_tsv(self.root / "TxDb-transcripts.tsv", fields, data)

    def check(self, model):
        with contextlib.redirect_stdout(io.StringIO()):
            return check_exports(model, self.roles, self.root, self.out, self.gff)

    def basic(self, gene_extra="Name=publicG", tx_extra="Name=publicT"):
        return self.model([row("gene", "g.original", extra=gene_extra),
                           row("mRNA", "t.original", "g.original", extra=tx_extra)])

    def test_reported_ID_Name_difference_and_three_TSS_examples(self):
        # IDs, Names and TSS below come from the user's three displayed rows.
        # End coordinates are synthetic; no full user annotation is available.
        g1, g2 = "Cre01.g000050_4532", "Cre01.g000100_4532"
        model = self.model([
            row("gene", g1 + ".v6.1", start=27912, end=30000, seqid="chromosome_01", extra="Name=" + g1),
            row("gene", g2 + ".v6.1", start=34228, end=50000, seqid="chromosome_01", extra="Name=" + g2),
            row("mRNA", g1 + ".1.v6.1", g1 + ".v6.1", 27912, 30000, seqid="chromosome_01", extra="Name=" + g1 + ".1"),
            row("mRNA", g2 + ".1.v6.1", g2 + ".v6.1", 34228, 40000, seqid="chromosome_01", extra="Name=" + g2 + ".1"),
            row("mRNA", g2 + ".2.v6.1", g2 + ".v6.1", 39579, 48000, seqid="chromosome_01", extra="Name=" + g2 + ".2")])
        self.observed([(g1 + ".1", "chromosome_01", "+", 27912, g1),
                       (g2 + ".1", "chromosome_01", "+", 34228, g2),
                       (g2 + ".2", "chromosome_01", "+", 39579, g2)])
        summary = self.check(model)
        self.assertEqual(summary["status"], "validated")
        self.assertEqual(summary["matched_associations"], 3)
        self.assertEqual(summary["transcripts_with_different_public_name"], 3)
        self.assertEqual(summary["genes_with_different_public_name"], 2)

    def test_unrelated_alias_strings_are_supported_without_species_rules(self):
        model = self.basic()
        self.observed([("publicT", "chr1", "+", 1, "publicG")])
        self.assertEqual(self.check(model)["matched_associations"], 1)
        mapping = (self.out / "TxDb-identity-map.tsv").read_text()
        self.assertIn("Name", mapping)

    def test_810_TE_aliases_preserve_original_relationships(self):
        records, observed = [], []
        for i in range(810):
            records += [row("transposable_element_gene", f"TE{i}", extra=f"Name=public-gene-{i}"),
                        row("mRNA", f"t{i}", f"TE{i}", extra=f"Name=RNA:{i}")]
            observed.append((f"RNA:{i}", "chr1", "+", 1, f"public-gene-{i}"))
        model = self.model(records)
        self.observed(observed)
        self.assertEqual(self.check(model)["matched_associations"], 810)

    def test_ID_fallback_without_names(self):
        model = self.basic("", "")
        self.observed([("t.original", "chr1", "+", 1, "g.original")])
        self.assertEqual(self.check(model)["transcripts_with_different_public_name"], 0)

    def test_explicit_gene_id_and_transcript_id(self):
        model = self.basic("Name=gene-name;gene_id=exportGene", "Name=tx-name;transcript_id=exportTx")
        self.observed([("exportTx", "chr1", "+", 1, "exportGene")])
        self.assertEqual(self.check(model)["status"], "validated")

    def test_suffix_is_never_removed_to_guess_alias(self):
        model = self.basic("", "")
        self.observed([("t", "chr1", "+", 1, "g")])
        with self.assertRaisesRegex(GFFError, "alias unresolved"):
            self.check(model)

    def test_duplicate_Name_fails_even_with_distinct_coordinates(self):
        model = self.model([row("gene", "g", end=500), row("mRNA", "a", "g", extra="Name=same"),
                            row("mRNA", "b", "g", 200, 300, extra="Name=same")])
        self.observed([("same", "chr1", "+", 1, "g"), ("same", "chr1", "+", 200, "g")])
        with self.assertRaisesRegex(GFFError, "alias unresolved"):
            self.check(model)

    def test_duplicate_gene_Name_cannot_merge_different_genes(self):
        model = self.model([row("gene", "g", extra="Name=same"), row("gene", "h", extra="Name=same"),
                            row("mRNA", "a", "g"), row("mRNA", "b", "h")])
        self.observed([("a", "chr1", "+", 1, "same"), ("b", "chr1", "+", 1, "same")])
        with self.assertRaisesRegex(GFFError, "gene alias unresolved"):
            self.check(model)

    def test_missing_transcript_and_extra_export_are_not_ignored(self):
        model = self.basic()
        self.observed([])
        with self.assertRaisesRegex(GFFError, "transcript missing"):
            self.check(model)
        self.observed([("publicT", "chr1", "+", 1, "publicG"), ("other", "chr1", "+", 1, "publicG")])
        with self.assertRaisesRegex(GFFError, "alias unresolved"):
            self.check(model)

    def test_wrong_TSS_sequence_or_strand_still_fails(self):
        model = self.basic()
        for seqid, strand, tss in (("chr1", "+", 2), ("chr2", "+", 1), ("chr1", "-", 1)):
            with self.subTest(seqid=seqid, strand=strand, tss=tss):
                self.observed([("publicT", seqid, strand, tss, "publicG")])
                with self.assertRaisesRegex(GFFError, "sequence/strand/TSS/span changed"):
                    self.check(model)

    def test_gene_reassignment_still_fails(self):
        model = self.model([row("gene", "g", extra="Name=G"), row("gene", "h", extra="Name=H"),
                            row("mRNA", "t", "g", extra="Name=T")])
        self.observed([("T", "chr1", "+", 1, "H")])
        with self.assertRaisesRegex(GFFError, "gene assignment changed"):
            self.check(model)

    def test_multi_parent_exact_set_required(self):
        model = self.model([row("gene", "g", extra="Name=G"), row("gene", "h", extra="Name=H"),
                            row("mRNA", "t", "g,h", extra="Name=T")])
        self.observed([("T", "chr1", "+", 1, "H,G")])
        self.assertEqual(self.check(model)["matched_associations"], 2)
        self.observed([("T", "chr1", "+", 1, "G")])
        with self.assertRaisesRegex(GFFError, "gene assignment changed"):
            self.check(model)

    def test_duplicate_TxDb_rows_are_not_accepted(self):
        model = self.basic()
        self.observed([("publicT", "chr1", "+", 1, "publicG")] * 2)
        with self.assertRaisesRegex(GFFError, "duplicate TxDb transcript"):
            self.check(model)

    def test_comma_in_gene_name_requires_lossless_edge_table(self):
        model = self.basic("Name=gene%2Cname", "Name=publicT")
        self.observed([("publicT", "chr1", "+", 1, "gene,name")])
        with self.assertRaisesRegex(GFFError, "legacy comma-joined"):
            self.check(model)
        self.observed([("publicT", "chr1", "+", 1, "gene,name")], edges=[(1, "gene,name")])
        self.assertEqual(self.check(model)["status"], "validated")

    def test_new_export_checks_full_transcript_span(self):
        model = self.basic()
        self.observed([("publicT", "chr1", "+", 1, "publicG")], edges=[(1, "publicG")])
        path = self.root / "TxDb-transcripts.tsv"
        with path.open(newline="") as handle:
            data = list(csv.DictReader(handle, delimiter="\t"))
        data[0].update(tx_start=1, tx_end=99)
        write_tsv(path, list(data[0]), data)
        with self.assertRaisesRegex(GFFError, "sequence/strand/TSS/span changed"):
            self.check(model)


if __name__ == "__main__":
    unittest.main()
