#!/usr/bin/env python3
"""Concordance check between synteny outputs produced at different miniprot
-M k-mer sampling exponents (see modules/local/miniprot_align.nf).

Reads the pipeline's own slider_links.tsv format directly (same file
bin/chain_blocks.mjs writes and plot_synteny_interactive.py reads --
target_chrom, target_start, target_end, reference_chrom, reference_start,
reference_end, score, orientation, mean_identity, anchor_density), so this
runs unmodified against real pipeline output, not a bespoke re-derivation.

One file is the baseline (normally the M=default run). Every other file is
scored against it: blocks are matched by (target_chrom, reference_chrom) plus a
reciprocal target-span overlap fraction -- deliberately NOT requiring an exact
coordinate match, since a different -M value chains a slightly different
anchor set and will not reproduce byte-identical block boundaries even when
it has found "the same" underlying synteny signal.

Usage:
    compare_links.py --baseline M1=path/to/M1.slider_links.tsv \\
                      --other    M2=path/to/M2.slider_links.tsv \\
                      --other    M3=path/to/M3.slider_links.tsv \\
                      [--overlap 0.5] [--out sweep_concordance.tsv]

Works identically on homeolog links (slider_homeolog_links.tsv) -- same
column schema, just a self-comparison.
"""
import argparse
import sys
from collections import defaultdict


def read_links(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = []
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            row = dict(zip(header, fields))
            row["target_start"] = int(row["target_start"])
            row["target_end"] = int(row["target_end"])
            row["reference_start"] = int(row["reference_start"])
            row["reference_end"] = int(row["reference_end"])
            row["score"] = int(row["score"])
            row["mean_identity"] = float(row["mean_identity"])
            row["anchor_density"] = float(row["anchor_density"])
            rows.append(row)
    return rows


def overlap_frac(a, b, start_key, end_key):
    lo = max(a[start_key], b[start_key])
    hi = min(a[end_key], b[end_key])
    if hi <= lo:
        return 0.0
    overlap = hi - lo
    span_a = a[end_key] - a[start_key]
    span_b = b[end_key] - b[start_key]
    return overlap / min(span_a, span_b) if min(span_a, span_b) > 0 else 0.0


def blocks_match(a, b, min_overlap):
    if a["target_chrom"] != b["target_chrom"] or a["reference_chrom"] != b["reference_chrom"]:
        return False
    q = overlap_frac(a, b, "target_start", "target_end")
    s = overlap_frac(a, b, "reference_start", "reference_end")
    return q >= min_overlap and s >= min_overlap


def match_greedy(baseline, other, min_overlap):
    """Greedy 1:1 match, largest-baseline-block-first. Returns
    (matched_pairs, baseline_only, other_only)."""
    by_pair = defaultdict(list)
    for j, b in enumerate(other):
        by_pair[(b["target_chrom"], b["reference_chrom"])].append(j)

    used_other = set()
    matched = []
    baseline_only = []
    for a in sorted(baseline, key=lambda r: -r["score"]):
        best_j, best_ov = None, 0.0
        for j in by_pair.get((a["target_chrom"], a["reference_chrom"]), []):
            if j in used_other:
                continue
            b = other[j]
            if not blocks_match(a, b, min_overlap):
                continue
            ov = min(overlap_frac(a, b, "target_start", "target_end"),
                     overlap_frac(a, b, "reference_start", "reference_end"))
            if ov > best_ov:
                best_j, best_ov = j, ov
        if best_j is not None:
            used_other.add(best_j)
            matched.append((a, other[best_j], best_ov))
        else:
            baseline_only.append(a)

    other_only = [other[j] for j in range(len(other)) if j not in used_other]
    return matched, baseline_only, other_only


def summarize(label, rows):
    if not rows:
        return dict(label=label, blocks=0, total_score=0, mean_identity=0.0)
    return dict(
        label=label,
        blocks=len(rows),
        total_score=sum(r["score"] for r in rows),
        mean_identity=sum(r["mean_identity"] * r["score"] for r in rows) / sum(r["score"] for r in rows),
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--baseline", required=True, metavar="LABEL=PATH",
                     help="the reference run, e.g. M1=path/to/links.tsv")
    ap.add_argument("--other", action="append", required=True, metavar="LABEL=PATH",
                     help="a run to score against the baseline; repeatable")
    ap.add_argument("--overlap", type=float, default=0.5,
                     help="min reciprocal query+subject overlap fraction to call two blocks "
                          "the same synteny call [0.5]")
    ap.add_argument("--out", default=None, help="also write the summary table here as TSV")
    args = ap.parse_args()

    def split(spec):
        label, _, path = spec.partition("=")
        if not path:
            sys.exit(f"ERROR: expected LABEL=PATH, got '{spec}'")
        return label, path

    base_label, base_path = split(args.baseline)
    baseline = read_links(base_path)
    base_summary = summarize(base_label, baseline)

    out_rows = []
    print(f"baseline: {base_label}  ({base_path})")
    print(f"  {base_summary['blocks']} blocks, {base_summary['total_score']} anchors, "
          f"mean_identity={base_summary['mean_identity']:.4f}")
    print()

    header = ["label", "blocks", "total_score", "mean_identity",
              "matched_vs_baseline", "match_rate", "baseline_only", "other_only",
              "recovered_score_frac"]
    print("\t".join(header))
    print("\t".join([base_label, str(base_summary["blocks"]), str(base_summary["total_score"]),
                      f"{base_summary['mean_identity']:.4f}", "-", "-", "-", "-", "-"]))
    out_rows.append([base_label, base_summary["blocks"], base_summary["total_score"],
                      f"{base_summary['mean_identity']:.4f}", "", "", "", "", ""])

    for spec in args.other:
        label, path = split(spec)
        other = read_links(path)
        summary = summarize(label, other)
        matched, baseline_only, other_only = match_greedy(baseline, other, args.overlap)

        match_rate = len(matched) / len(baseline) if baseline else 0.0
        recovered_score = sum(a["score"] for a, _, _ in matched)
        recovered_frac = recovered_score / base_summary["total_score"] if base_summary["total_score"] else 0.0

        row = [label, summary["blocks"], summary["total_score"], f"{summary['mean_identity']:.4f}",
               len(matched), f"{match_rate:.3f}", len(baseline_only), len(other_only),
               f"{recovered_frac:.3f}"]
        print("\t".join(str(x) for x in row))
        out_rows.append(row)

    if args.out:
        with open(args.out, "w") as f:
            f.write("\t".join(header) + "\n")
            for row in out_rows:
                f.write("\t".join(str(x) for x in row) + "\n")
        print(f"\nwritten: {args.out}", file=sys.stderr)

    print(
        "\nReading this: match_rate = fraction of the BASELINE's blocks that "
        "still show up (by overlap) at this -M value -- the number to quote "
        "as \"% of synteny calls retained.\" recovered_score_frac is the same "
        "idea weighted by block support (anchor count), so a few small lost "
        "blocks don't look as bad as losing a few huge, well-supported ones.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
