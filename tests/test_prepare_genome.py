#!/usr/bin/env python3
"""Equivalence test for the one-pass rewrite of bin/rename_sequences.py
(plan A: pipeline I/O).

Compares the new script's output against the two scripts it replaces --
the old bin/rename_sequences.py (rename + write a renamed FASTA) and the old
bin/find_assembly_gaps.py (gap scan over that renamed FASTA), both fetched
straight from git at commit eac2c86 so there's no risk of testing against a
stale local copy -- plus faSize -detailed, run through the ucsc-fasize
container per docs/plans/README.md, as the oracle for --out_sizes.

Plain asserts, no pytest, no third-party libraries: this has to run
unmodified under both the host's Python 3.10 and the pipeline containers'
Python 3.13. Needs `git` and a working `docker` (for the faSize oracle only)
on PATH.
"""
import gzip
import os
import subprocess
import sys
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIXTURES_SCRIPT = os.path.join(REPO_ROOT, 'tests', 'fixtures', 'make_mini_genome.py')
NEW_RENAME_PY = os.path.join(REPO_ROOT, 'bin', 'rename_sequences.py')
OLD_COMMIT = 'eac2c86'
FASIZE_IMAGE = 'quay.io/biocontainers/ucsc-fasize:482--h0b57e2e_0'
SEED = 20260923
VARIANTS = ['scan.fa', 'scan_crlf.fa', 'scan.fa.gz']
MIN_GAPS = [1, 10, 100]
MIN_SEQ_SIZES = [0, 1000]


def run(cmd, **kwargs):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kwargs)


def git_show(relpath, dest_path):
    """Materialize a file exactly as it was at OLD_COMMIT, as the oracle --
    not whatever bin/ currently holds after this plan's edits."""
    result = run(['git', 'show', f'{OLD_COMMIT}:{relpath}'], cwd=REPO_ROOT)
    with open(dest_path, 'w') as f:
        f.write(result.stdout)
    os.chmod(dest_path, 0o755)
    return dest_path


def read(path):
    with open(path) as f:
        return f.read()


def run_old_rename(old_rename_py, fasta_path, workdir, tag):
    out_fasta = os.path.join(workdir, f'{tag}.old.renamed.fa')
    out_lookup = os.path.join(workdir, f'{tag}.old.lookup.tsv')
    run([sys.executable, old_rename_py,
         '--fasta', fasta_path, '--out_fasta', out_fasta, '--out_lookup', out_lookup])
    return out_fasta, out_lookup


def run_old_gaps(old_gaps_py, renamed_fasta, min_gap, workdir, tag):
    out_gaps = os.path.join(workdir, f'{tag}.old.gaps.{min_gap}.tsv')
    run([sys.executable, old_gaps_py,
         '--fasta', renamed_fasta, '--min_gap', str(min_gap), '--out', out_gaps])
    return out_gaps


def run_new_rename(fasta_path, min_seq_size, min_gap, workdir, tag):
    out_lookup = os.path.join(workdir, f'{tag}.new.lookup.{min_seq_size}.{min_gap}.tsv')
    out_sizes = os.path.join(workdir, f'{tag}.new.sizes.{min_seq_size}.{min_gap}.tsv')
    out_gaps = os.path.join(workdir, f'{tag}.new.gaps.{min_seq_size}.{min_gap}.tsv')
    run([sys.executable, NEW_RENAME_PY,
         '--fasta', fasta_path,
         '--out_lookup', out_lookup,
         '--out_sizes', out_sizes, '--min_seq_size', str(min_seq_size),
         '--out_gaps', out_gaps, '--min_gap', str(min_gap)])
    return out_lookup, out_sizes, out_gaps


def fasize_detailed(renamed_fasta):
    """faSize -detailed <fasta>: 'seqid\\tlength' per line, FASTA order, no
    summary line (checked directly against the container -- see the plan's
    acceptance check 1)."""
    d = os.path.dirname(os.path.abspath(renamed_fasta))
    fname = os.path.basename(renamed_fasta)
    result = run(['docker', 'run', '--rm', '--platform', 'linux/amd64',
                  '-v', f'{d}:{d}', '-w', d, FASIZE_IMAGE,
                  'faSize', '-detailed', fname])
    sizes = []
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        seq_id, length = line.split('\t')
        sizes.append((seq_id, int(length)))
    return sizes


def main():
    with tempfile.TemporaryDirectory(prefix='test_prepare_genome_') as workdir:
        print(f"[test] workdir: {workdir}")

        old_rename_py = git_show('bin/rename_sequences.py', os.path.join(workdir, 'old_rename_sequences.py'))
        old_gaps_py = git_show('bin/find_assembly_gaps.py', os.path.join(workdir, 'old_find_assembly_gaps.py'))

        fixtures_dir = os.path.join(workdir, 'fixtures')
        os.makedirs(fixtures_dir)
        run([sys.executable, FIXTURES_SCRIPT, '--seed', str(SEED), '--outdir', fixtures_dir])

        n_checks = 0
        for variant in VARIANTS:
            variant_path = os.path.join(fixtures_dir, variant)

            if variant.endswith('.gz'):
                # the old script can't read gzip -- compare against its
                # output on the decompressed file (see the plan)
                old_input_path = os.path.join(workdir, f'{variant}.decompressed.fa')
                with gzip.open(variant_path, 'rb') as fin, open(old_input_path, 'wb') as fout:
                    fout.write(fin.read())
            else:
                old_input_path = variant_path

            old_renamed_fa, old_lookup_path = run_old_rename(old_rename_py, old_input_path, workdir, variant)
            old_lookup = read(old_lookup_path)

            all_sizes = fasize_detailed(old_renamed_fa)

            for min_gap in MIN_GAPS:
                old_gaps_path = run_old_gaps(old_gaps_py, old_renamed_fa, min_gap, workdir, variant)
                old_gaps = read(old_gaps_path)

                for min_seq_size in MIN_SEQ_SIZES:
                    new_lookup_path, new_sizes_path, new_gaps_path = run_new_rename(
                        variant_path, min_seq_size, min_gap, workdir, variant)

                    ctx = f"{variant} min_gap={min_gap} min_seq_size={min_seq_size}"

                    new_lookup = read(new_lookup_path)
                    assert new_lookup == old_lookup, f"lookup mismatch ({ctx})"
                    n_checks += 1

                    new_gaps = read(new_gaps_path)
                    assert new_gaps == old_gaps, f"gaps mismatch ({ctx})"
                    n_checks += 1

                    expected_sizes = ''.join(f"{seq_id}\t{length}\n"
                                              for seq_id, length in all_sizes
                                              if length >= min_seq_size)
                    new_sizes = read(new_sizes_path)
                    assert new_sizes == expected_sizes, f"sizes mismatch vs faSize|awk ({ctx})"
                    n_checks += 1

            print(f"[test] {variant}: lookup/gaps/sizes match across "
                  f"{len(MIN_GAPS)}x{len(MIN_SEQ_SIZES)} (min_gap, min_seq_size) combinations")

        print(f"[test] ALL {n_checks} CHECKS PASSED "
              f"({len(VARIANTS)} input variants x {len(MIN_GAPS)} min_gap x {len(MIN_SEQ_SIZES)} min_seq_size, "
              f"x{{lookup,gaps,sizes}})")


if __name__ == '__main__':
    main()
