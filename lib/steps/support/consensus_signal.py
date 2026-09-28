#!/usr/bin/env python3
"""Fill consensus BED scores with the maximum overlapping MACS signalValue.
nf-core consensus BED column 5 is a placeholder. Source values come from
narrowPeak column 7, combined across replicates. No bedtools dependency.
Usage: consensus_signal.py <consensus.bed> <all_narrowPeak.bed> <out.bed>
Output: chr start end name signalValue strand
"""

import bisect
import collections
import sys


def load_narrowpeak(path):
    """Collect (start, end, signalValue) per chromosome, sorted by start."""
    by_chr = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 7:
                continue
            try:
                st, en, sig = int(float(c[1])), int(float(c[2])), float(c[6])
            except ValueError:
                continue          # 表头或坏行
            by_chr[c[0]].append((st, en, sig))
    for ch in by_chr:
        by_chr[ch].sort()
    return by_chr


def build_index(intervals):
    """Index starts and prefix maximum ends, without a maximum peak length."""
    starts = [x[0] for x in intervals]
    prefix_ends = []
    for _, end, _ in intervals:
        prefix_ends.append(max(end, prefix_ends[-1] if prefix_ends else end))
    return starts, prefix_ends, intervals


def max_signal(starts, prefix_ends, intervals, qs, qe):
    """Return the maximum signal overlapping [qs, qe), or None."""
    # 所有可能重叠的区间，其 start < qe
    hi = bisect.bisect_left(starts, qe)
    best = None
    for i in range(hi - 1, -1, -1):
        if prefix_ends[i] <= qs:
            break
        st, en, sig = intervals[i]
        if en <= qs:
            continue
        if best is None or sig > best:
            best = sig
    return best


def main():
    if len(sys.argv) < 4:
        sys.exit("Usage: consensus_signal.py <consensus.bed> <all_narrowPeak.bed> <out.bed>")
    cons_f, np_f, out_f = sys.argv[1], sys.argv[2], sys.argv[3]

    by_chr = load_narrowpeak(np_f)
    if not by_chr:
        sys.exit("[ERROR] narrowPeak has no usable intervals")
    index = {ch: build_index(v) for ch, v in by_chr.items()}

    n_out = n_hit = 0
    rows = []
    with open(cons_f) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 3:
                continue
            try:
                qs, qe = int(float(c[1])), int(float(c[2]))
            except ValueError:
                continue          # 表头
            ch = c[0]
            name = c[3] if len(c) > 3 and c[3] else "peak"
            sig = None
            if ch in index:
                sig = max_signal(*index[ch], qs, qe)
            if sig is not None:
                n_hit += 1
            rows.append((ch, qs, qe, name, sig if sig is not None else 0.0))
            n_out += 1

    if n_out == 0:
        sys.exit("[ERROR] Consensus file has no usable intervals")

    rows.sort(key=lambda r: (r[0], r[1]))
    with open(out_f, "w") as out:
        for ch, qs, qe, name, sig in rows:
            out.write("%s\t%d\t%d\t%s\t%g\t+\n" % (ch, qs, qe, name, sig))

    print("  Consensus signalValue: %d intervals, %d matched to a score" % (n_out, n_hit),
          file=sys.stderr)
    if n_hit == 0:
        sys.exit("[ERROR] No consensus intervals match narrowPeak; check sequence names")


if __name__ == "__main__":
    main()
