#!/usr/bin/env python3
"""
Cell-type-aware UMI annotation of SComatic Step 2 PASS variants.

For each Step 2 record, examine its supporting cell types using their
respective BAMs. Preserve original SComatic cell-type ordering.

Outputs:
  1. Cell-type-specific UMI summary
  2. SComatic-compatible variant summary (comma-aligned UMI fields)
  3. Per-molecule read evidence

Molecules are defined by (CB, UB) within each variant and cell type.
Only reads with valid corrected UB tags contribute to UMI counts.
"""

import argparse
import csv
import gzip
from collections import Counter, defaultdict
from contextlib import ExitStack
from pathlib import Path

import pysam


UMI_FIELDS = [
    "REF_CONSENSUS_UMIS",
    "ALT_CONSENSUS_UMIS",
    "UMI_VAF_STRICT",
    "ALT_SUPPORTING_CELLS",
    "CONFLICTING_UMIS",
    "OTHER_OR_MIXED_UMIS",
    "READS_MISSING_UB",
    "READS_MISSING_CB",
    "POSSIBLY_DEPTH_CAPPED",
]


UMI_OUTPUT_NAMES = {
    "REF_CONSENSUS_UMIS": "UMIs_WITH_ONLY_REF_READS",
    "ALT_CONSENSUS_UMIS": "UMIs_WITH_ONLY_ALT_READS",
    "UMI_VAF_STRICT": "ALT_FRACTION_OF_REF_ALT_CONSISTENT_UMIs",
    "ALT_SUPPORTING_CELLS": "CELLS_WITH_ALT_SUPPORTING_UMIs",
    "CONFLICTING_UMIS": "UMIs_WITH_BOTH_REF_AND_ALT_READS",
    "OTHER_OR_MIXED_UMIS": "UMIs_WITH_OTHER_OR_DISCORDANT_READS",
    "READS_MISSING_UB": "READS_WITHOUT_CORRECTED_UMI_TAG",
    "READS_MISSING_CB": "READS_WITHOUT_CELL_BARCODE",
    "POSSIBLY_DEPTH_CAPPED": "PILEUP_REACHED_MAX_DEPTH",
}

CELL_FIELDS = [
    "SAMPLE", "RECORD_INDEX", "CELL_TYPE", "CELL_TYPE_INDEX",
    "CHROM", "POS", "REF", "ALT", "VARIANT_KEY",
    "SCOMATIC_Bc", "SCOMATIC_Cc", "SCOMATIC_Dp",
    "SCOMATIC_VAF", "SCOMATIC_CELL_TYPE_FILTER",
] + [UMI_OUTPUT_NAMES[field] for field in UMI_FIELDS]

MOLECULE_FIELDS = [
    "SAMPLE", "RECORD_INDEX", "CELL_TYPE", "VARIANT_KEY",
    "CB", "UB", "N_READS",
    "REF_READS", "ALT_READS", "OTHER_READS", "STATUS",
]


def read_step2(path, max_records):
    """Read Step 2 records, ignoring SComatic metadata lines."""
    opener = gzip.open if str(path).endswith(".gz") else open

    with opener(path, "rt", newline="") as handle:
        lines = (
            line for line in handle
            if line.strip() and not line.startswith("##")
        )

        reader = csv.DictReader(lines, delimiter="\t")

        if not reader.fieldnames:
            raise ValueError("No Step 2 header found")

        reader.fieldnames = [
            name.lstrip("#") for name in reader.fieldnames
        ]

        required = {
            "CHROM", "Start", "REF", "ALT",
            "FILTER", "Cell_types"
        }
        missing = required - set(reader.fieldnames)

        if missing:
            raise ValueError(f"Missing Step 2 columns: {missing}")

        columns = list(reader.fieldnames)
        records = []

        for row in reader:
            if row["FILTER"].strip() != "PASS":
                continue

            records.append(row)

            if max_records is not None and len(records) >= max_records:
                break

    return columns, records


def parse_cell_types(row):
    """
    Parse aligned Step 2 cell-type-specific fields.

    ALT, Cell_types and the SComatic per-cell-type measurements
    are positionally aligned.
    """
    cell_types = [
        x.strip() for x in row["Cell_types"].split(",")
    ]

    if not cell_types or any(not x for x in cell_types):
        raise ValueError(f"Invalid Cell_types: {row['Cell_types']}")

    alts = [x.strip().upper() for x in row["ALT"].split(",")]

    if len(alts) == 1:
        alts = alts * len(cell_types)
    elif len(alts) != len(cell_types):
        raise ValueError(
            f"ALT/Cell_types mismatch: "
            f"{row['ALT']} / {row['Cell_types']}"
        )

    def aligned(field):
        text = row.get(field, "") or ""

        if not text:
            return [""] * len(cell_types)

        values = [x.strip() for x in text.split(",")]

        if len(values) == 1:
            return values * len(cell_types)

        if len(values) != len(cell_types):
            raise ValueError(
                f"{field} has {len(values)} values but "
                f"Cell_types has {len(cell_types)}: {text}"
            )

        return values

    bc = aligned("Bc")
    cc = aligned("Cc")
    dp = aligned("Dp")
    vaf = aligned("VAF")
    filters = aligned("Cell_type_Filter")

    return [
        {
            "cell_type": cell_type,
            "cell_type_index": i,
            "alt": alts[i],
            "bc": bc[i],
            "cc": cc[i],
            "dp": dp[i],
            "vaf": vaf[i],
            "filter": filters[i],
        }
        for i, cell_type in enumerate(cell_types)
    ]


def analyse_site(
    bam, chrom, pos, ref, alt, args,
    sample, cell_type, record_index, molecule_writer
):
    """Calculate read and UMI evidence at one SNV in one BAM."""

    if chrom not in bam.references:
        raise RuntimeError(
            f"{chrom} is absent from BAM for {cell_type}"
        )

    molecules = defaultdict(Counter)
    alleles = Counter()

    read_depth = 0
    umi_read_depth = 0
    missing_cb = 0
    missing_ub = 0
    possibly_capped = False

    pileup = bam.pileup(
        chrom,
        pos - 1,
        pos,
        truncate=True,
        stepper="nofilter",
        min_base_quality=0,
        min_mapping_quality=0,
        max_depth=args.max_depth,
        ignore_overlaps=False,
    )

    for col in pileup:
        if col.reference_pos != pos - 1:
            continue

        if col.nsegments >= args.max_depth:
            possibly_capped = True

        for pr in col.pileups:
            read = pr.alignment

            if (
                read.is_unmapped
                or read.is_secondary
                or read.is_supplementary
                or read.is_duplicate
                or read.is_qcfail
                or read.mapping_quality < args.min_mq
                or pr.is_del
                or pr.is_refskip
                or pr.query_position is None
            ):
                continue

            q = pr.query_position

            if (
                read.query_qualities is None
                or read.query_sequence is None
                or read.query_qualities[q] < args.min_bq
            ):
                continue

            base = read.query_sequence[q].upper()

            if base not in "ACGT":
                continue

            if not read.has_tag("CB"):
                missing_cb += 1
                continue

            cb = str(read.get_tag("CB")).split("-")[0]

            if not cb:
                missing_cb += 1
                continue

            read_depth += 1
            alleles[base] += 1

            if not read.has_tag("UB"):
                missing_ub += 1
                continue

            ub = str(read.get_tag("UB"))

            if not ub or ub in (".", "NA", "None"):
                missing_ub += 1
                continue

            umi_read_depth += 1
            molecules[(cb, ub)][base] += 1

    ref_umis = 0
    alt_umis = 0
    ref_consensus = 0
    alt_consensus = 0
    single_read = 0
    multiread = 0
    multiread_alt = 0
    multiread_ref = 0
    conflicting = 0
    mixed = 0

    alt_cells = set()
    variant_key = f"{chrom}:{pos}:{ref}:{alt}"

    for (cb, ub), counts in molecules.items():
        r = counts[ref]
        a = counts[alt]
        o = sum(
            n for base, n in counts.items()
            if base not in (ref, alt)
        )
        n = sum(counts.values())

        ref_umis += int(r > 0)
        alt_umis += int(a > 0)

        single_read += int(n == 1)
        multiread += int(n >= 2)

        if a:
            alt_cells.add(cb)

        if r and a:
            conflicting += 1

        if a and not r and not o:
            status = "CONSENSUS_ALT"
            alt_consensus += 1

            if n >= 2:
                multiread_alt += 1

        elif r and not a and not o:
            status = "CONSENSUS_REF"
            ref_consensus += 1

            if n >= 2:
                multiread_ref += 1

        else:
            status = "MIXED_OR_OTHER"
            mixed += 1

        molecule_writer.writerow({
            "SAMPLE": sample,
            "RECORD_INDEX": record_index,
            "CELL_TYPE": cell_type,
            "VARIANT_KEY": variant_key,
            "CB": cb,
            "UB": ub,
            "N_READS": n,
            "REF_READS": r,
            "ALT_READS": a,
            "OTHER_READS": o,
            "STATUS": status,
        })

    denominator = ref_consensus + alt_consensus

    return {
        "READ_DEPTH": read_depth,
        "UMI_READ_DEPTH": umi_read_depth,
        "UMI_DEPTH": len(molecules),
        "REF_READS": alleles[ref],
        "ALT_READS": alleles[alt],
        "REF_UNIQUE_UMIS": ref_umis,
        "ALT_UNIQUE_UMIS": alt_umis,
        "REF_CONSENSUS_UMIS": ref_consensus,
        "ALT_CONSENSUS_UMIS": alt_consensus,
        "UMI_VAF_STRICT": (
            f"{alt_consensus / denominator:.6f}"
            if denominator else "NA"
        ),
        "ALT_SUPPORTING_CELLS": len(alt_cells),
        "SINGLE_READ_UMIS": single_read,
        "MULTIREAD_UMIS": multiread,
        "MULTIREAD_ALT_CONSENSUS_UMIS": multiread_alt,
        "MULTIREAD_REF_CONSENSUS_UMIS": multiread_ref,
        "CONFLICTING_UMIS": conflicting,
        "OTHER_OR_MIXED_UMIS": mixed,
        "READS_MISSING_UB": missing_ub,
        "READS_MISSING_CB": missing_cb,
        "POSSIBLY_DEPTH_CAPPED": possibly_capped,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--step2", required=True)
    parser.add_argument("--bam-dir", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out-prefix", required=True)

    parser.add_argument("--max-records", type=int, default=None)
    parser.add_argument("--min-bq", type=int, default=20)
    parser.add_argument("--min-mq", type=int, default=255)
    parser.add_argument("--max-depth", type=int, default=8000)

    args = parser.parse_args()

    if args.max_records is not None and args.max_records < 1:
        parser.error("--max-records must be positive")

    if args.max_depth < 1:
        parser.error("--max-depth must be positive")

    bam_dir = Path(args.bam_dir)
    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    cell_output = Path(str(prefix) + ".umi_celltype_summary.tsv")
    variant_output = Path(str(prefix) + ".umi_variant_summary.tsv")
    molecule_output = Path(str(prefix) + ".umi_molecules.tsv")

    original_columns, records = read_step2(
        args.step2, args.max_records
    )

    if not records:
        raise RuntimeError("No Step 2 PASS records found")

    print(
        f"Selected {len(records)} Step 2 PASS records "
        f"for {args.sample}",
        flush=True,
    )

    # Preserve every original Step 2 column.
    # New UMI fields use a prefix to avoid overwriting SComatic.
    variant_columns = (
        ["#CHROM" if col == "CHROM" else col for col in original_columns]
        + ["SAMPLE", "RECORD_INDEX"]
        + [UMI_OUTPUT_NAMES[field] for field in UMI_FIELDS]
    )

    bam_handles = {}

    with ExitStack() as stack:
        cell_handle = stack.enter_context(
            open(cell_output, "w", newline="")
        )
        variant_handle = stack.enter_context(
            open(variant_output, "w", newline="")
        )
        molecule_handle = stack.enter_context(
            open(molecule_output, "w", newline="")
        )

        cell_writer = csv.DictWriter(
            cell_handle, fieldnames=CELL_FIELDS, delimiter="\t"
        )
        variant_writer = csv.DictWriter(
            variant_handle, fieldnames=variant_columns,
            delimiter="\t"
        )
        molecule_writer = csv.DictWriter(
            molecule_handle, fieldnames=MOLECULE_FIELDS,
            delimiter="\t"
        )

        cell_writer.writeheader()
        variant_writer.writeheader()
        molecule_writer.writeheader()

        for record_index, row in enumerate(records, 1):
            chrom = row["CHROM"]
            pos = int(row["Start"])
            ref = row["REF"].upper()

            entries = parse_cell_types(row)

            print(
                f"[{record_index}/{len(records)}] "
                f"{chrom}:{pos}:{ref} "
                f"cell_types={row['Cell_types']}",
                flush=True,
            )

            # Each field is accumulated in the exact Cell_types order.
            umi_values = {field: [] for field in UMI_FIELDS}

            for entry in entries:
                cell_type = entry["cell_type"]
                alt = entry["alt"]

                if (
                    len(ref) != 1
                    or len(alt) != 1
                    or ref not in "ACGT"
                    or alt not in "ACGT"
                    or ref == alt
                ):
                    raise ValueError(
                        f"Invalid SNV at record {record_index}: "
                        f"{chrom}:{pos}:{ref}:{alt}"
                    )

                if cell_type not in bam_handles:
                    bam_path = (
                        bam_dir /
                        f"{args.sample}.{cell_type}.bam"
                    )

                    if not bam_path.is_file():
                        raise FileNotFoundError(
                            f"Missing BAM: {bam_path}"
                        )

                    print(
                        f"  Opening {bam_path.name}",
                        flush=True,
                    )

                    bam_handles[cell_type] = stack.enter_context(
                        pysam.AlignmentFile(str(bam_path), "rb")
                    )

                metrics = analyse_site(
                    bam=bam_handles[cell_type],
                    chrom=chrom,
                    pos=pos,
                    ref=ref,
                    alt=alt,
                    args=args,
                    sample=args.sample,
                    cell_type=cell_type,
                    record_index=record_index,
                    molecule_writer=molecule_writer,
                )

                cell_writer.writerow({
                    "SAMPLE": args.sample,
                    "RECORD_INDEX": record_index,
                    "CELL_TYPE": cell_type,
                    "CELL_TYPE_INDEX": entry["cell_type_index"],
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt,
                    "VARIANT_KEY": f"{chrom}:{pos}:{ref}:{alt}",
                    "SCOMATIC_Bc": entry["bc"],
                    "SCOMATIC_Cc": entry["cc"],
                    "SCOMATIC_Dp": entry["dp"],
                    "SCOMATIC_VAF": entry["vaf"],
                    "SCOMATIC_CELL_TYPE_FILTER": entry["filter"],
                    **{
                        UMI_OUTPUT_NAMES[field]: metrics[field]
                        for field in UMI_FIELDS
                    },
                })

                for field in UMI_FIELDS:
                    umi_values[field].append(str(metrics[field]))

                print(
                    f"  {cell_type}: "
                    f"reads={metrics['READ_DEPTH']} "
                    f"UMIs={metrics['UMI_DEPTH']} "
                    f"ALT_strict={metrics['ALT_CONSENSUS_UMIS']} "
                    f"VAF={metrics['UMI_VAF_STRICT']}",
                    flush=True,
                )

            variant_row = dict(row)
            variant_row["#CHROM"] = variant_row.pop("CHROM")
            variant_row["SAMPLE"] = args.sample
            variant_row["RECORD_INDEX"] = record_index

            for field in UMI_FIELDS:
                output_field = UMI_OUTPUT_NAMES[field]
                variant_row[output_field] = ",".join(umi_values[field])

            variant_writer.writerow(variant_row)

            cell_handle.flush()
            variant_handle.flush()
            molecule_handle.flush()

    print("\nCOMPLETE", flush=True)
    print(f"Cell-type summary: {cell_output}", flush=True)
    print(f"Variant summary:   {variant_output}", flush=True)
    print(f"Molecule evidence: {molecule_output}", flush=True)


if __name__ == "__main__":
    main()
