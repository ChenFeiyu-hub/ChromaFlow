# FIMO zero-hit compatibility hotfix for ChromaFlow 1.2.8

Patch identifier: `fimo-zero-hits-v1` (2026-11-06). The CLI version remains 1.2.8;
the installed manifest and batch `[HOTFIX]` message identify this patch.

FIMO 5.5.7 prints the TSV header only when it prints its first reported match.
If no matches pass the scan threshold, `fimo.tsv` contains only a blank line and
the three version/format/command comments. This is a valid zero-hit result.
The original ChromaFlow reader rejected it as missing required columns.
See the official 5.5.7 source archive's `src/fimo-output.c`,
`print_site_as_tsv()` and `print_fimo_tsv_file()`:
https://meme-suite.org/meme/meme-software/5.5.7/meme-5.5.7.tar.gz

The corrected reader accepts that exact completion-footer shape and records zero
hits for the motif. Other motifs remain eligible under the unchanged filters.
Blank files, incomplete/unrecognized comments, malformed tables, and contradictory
GFF hits remain errors. A validated companion GFF is included in selection provenance.
No FASTQ, FASTA, annotation, FIMO result, gene identifier or threshold is edited.

After the active analysis has ended, transfer the hotfix archive and checksum to
the server. Run as the owner of the ChromaFlow installation:

```bash
cd ~/nfcore-upload
sha256sum -c ChromaFlow-1.2.8-fimo-zero-hits-hotfix.tar.gz.sha256 &&
tar -xzf ChromaFlow-1.2.8-fimo-zero-hits-hotfix.tar.gz &&
python3 ChromaFlow-1.2.8-fimo-zero-hits-hotfix/install.py "$HOME/biosoft/ChromaFlow-1.2.8"
```

The installer checks old and new file hashes, backs up replaced files under the
installation's `.chromaflow-hotfixes/` directory, and refuses unexpected changes.
Its manifest records the entire runtime tree before and after the update.
No Conda environment rebuild or new container download is required.

Resume with the same project config, sample table and selection as the failed run:

```bash
CF="$HOME/biosoft/ChromaFlow-1.2.8/bin/chromaflow"
"$CF" batch run -c chromaflow.config --samples samples_dof_smoke.tsv --resume --until candidates
"$CF" batch status -c chromaflow.config
```

When inputs, environments and outputs still match, reference preparation, nf-core,
annotation and motif discovery are reused. Candidate integration is recomputed for
all selected analyses under the updated reader, including previously successful
ones. A subsequent unchanged resume reuses those updated candidate results too.
Look for `[RESUME]` upstream and `[FIMO] No reported hits (completed footer-only
result)` in the affected candidate log. `fimo_out_5` is a run-local directory index,
not a stable motif name. Use the raw ID or motif catalog for explicit selection.

Later runtime edits disable this narrow compatibility mapping; normal conservative
cache invalidation applies. Never edit cached success records to bypass checks.
The separate 130 bp nf-core parameter compatibility issue is not fixed by this patch.
