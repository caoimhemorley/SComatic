#!/usr/bin/env python3
"""Safely merge coordinate-sorted SComatic cell-type BAMs for one sample."""
import argparse
import glob
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def run(args):
    print('+', ' '.join(map(str, args)), flush=True)
    subprocess.run(args, check=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample', required=True)
    p.add_argument('--input-dir', required=True, help='Directory recursively containing split BAMs')
    p.add_argument('--output', required=True, help='Output merged BAM path')
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--expected-cell-types', type=int, help='Optional expected number of split BAMs')
    p.add_argument('--overwrite', action='store_true')
    p.add_argument('--check-overlap', action='store_true', help='Expensive full alignment-overlap check; may require substantial temporary disk space')
    args = p.parse_args()

    if not shutil.which('samtools'):
        p.error('samtools is required on PATH')
    if args.threads < 1:
        p.error('--threads must be >= 1')

    source = Path(args.input_dir).resolve()
    output = Path(args.output).resolve()
    if not source.is_dir():
        p.error(f'Input directory does not exist: {source}')
    pattern = f'{glob.escape(args.sample)}.*.bam'
    bams = sorted(x.resolve() for x in source.rglob(pattern) if x.is_file() and not x.name.endswith('.bam.bai'))
    if not bams:
        p.error(f'No split BAMs matching {pattern} found under {source}')
    if output in bams:
        p.error('Output path is one of the input BAMs')
    if output.exists() and not args.overwrite:
        p.error(f'Output exists: {output} (use --overwrite to replace)')
    if args.expected_cell_types is not None and len(bams) != args.expected_cell_types:
        p.error(f'Expected {args.expected_cell_types} BAMs but found {len(bams)}')

    labels = [x.name[len(args.sample) + 1:-4] for x in bams]
    if len(set(labels)) != len(labels):
        p.error('Duplicate cell-type labels detected; possible multiple job attempts or copies. Select a directory containing exactly one final BAM per cell type.')

    print(f'Sample: {args.sample}\nInput BAMs ({len(bams)}):', flush=True)
    for label, bam in zip(labels, bams):
        print(f'  {label}: {bam}', flush=True)

    try:
        import pysam
    except ImportError:
        p.error('pysam is required to verify headers and check for duplicate alignment records')

    reference = None
    read_groups = {}
    for bam in bams:
        with pysam.AlignmentFile(str(bam), 'rb') as f:
            hd = f.header.to_dict().get('HD', {})
            if hd.get('SO') != 'coordinate':
                p.error(f'Not declared coordinate-sorted: {bam}')
            sq = tuple((r['SN'], r['LN']) for r in f.header.to_dict().get('SQ', []))
            if reference is None:
                reference = sq
            elif sq != reference:
                p.error(f'Reference sequence dictionaries differ: {bam}')
            for rg in f.header.to_dict().get('RG', []):
                rgid = rg.get('ID')
                if rgid and rgid in read_groups and read_groups[rgid] != rg:
                    p.error(f'Conflicting read-group ID {rgid} between BAM headers')
                if rgid:
                    read_groups[rgid] = rg

    if args.check_overlap:
        # Verify disjointness of CB-tagged alignments without holding all read names in memory.
        # A read can legitimately have multiple alignments within one BAM, so compare across BAMs.
        # SQLite stores the keys on disk to keep RAM usage bounded.
        import sqlite3
        with tempfile.TemporaryDirectory(prefix='scomatic_merge_check_') as tmp:
            db = sqlite3.connect(str(Path(tmp) / 'read_keys.sqlite'))
            db.execute('PRAGMA journal_mode=OFF')
            db.execute('PRAGMA synchronous=OFF')
            db.execute('CREATE TABLE reads (key TEXT PRIMARY KEY, source INTEGER)')
            duplicate_examples = []
            total = 0
            for i, bam in enumerate(bams):
                with pysam.AlignmentFile(str(bam), 'rb') as f:
                    for read in f.fetch(until_eof=True):
                        # Only exact same alignment across separate BAMs is an error.
                        # Include flag, position, mate position, CIGAR, and CB.
                        key = '\x1f'.join((read.query_name or '', str(read.flag), str(read.reference_id),
                                           str(read.reference_start), str(read.next_reference_id),
                                           str(read.next_reference_start), read.cigarstring or '',
                                           read.get_tag('CB') if read.has_tag('CB') else ''))
                        row = db.execute('SELECT source FROM reads WHERE key=?', (key,)).fetchone()
                        if row is None:
                            db.execute('INSERT INTO reads VALUES (?, ?)', (key, i))
                        elif row[0] != i:
                            duplicate_examples.append((labels[row[0]], labels[i], read.query_name))
                            if len(duplicate_examples) >= 5:
                                break
                        total += 1
                db.commit()
                if duplicate_examples:
                    break
                print(f'Checked {labels[i]}: {i + 1}/{len(bams)} BAMs', flush=True)
            db.close()
            if duplicate_examples:
                p.error(f'Identical alignments found across cell-type BAMs (first examples: {duplicate_examples}). Aborting to avoid double-counting.')

    if not args.check_overlap:
        print('Note: full cross-BAM alignment-overlap check skipped; use --check-overlap for validation.', flush=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.merge_scomatic_', dir=str(output.parent)) as tmp:
        merged = Path(tmp) / 'merged.bam'
        run(['samtools', 'merge', '-@', str(args.threads), '-o', str(merged), *map(str, bams)])
        run(['samtools', 'quickcheck', '-v', str(merged)])
        run(['samtools', 'index', '-@', str(args.threads), str(merged)])
        if output.exists() and not args.overwrite:
            p.error(f'Output appeared during merge: {output}')
        os.replace(merged, output)
        os.replace(Path(str(merged) + '.bai'), Path(str(output) + '.bai'))
    print(f'Completed: {output}\nIndex: {output}.bai', flush=True)


if __name__ == '__main__':
    main()
