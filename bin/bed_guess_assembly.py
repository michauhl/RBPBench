#!/usr/bin/env python3

import argparse
import csv
import os
import subprocess
import tempfile
from pathlib import Path

from rbpbench import benchlib

"""
Guess assembly of a given BED regions file.

Example call:
python bed_guess_assembly.py --in regions_to_guess.bed --mode all --bed-sc-thr 2
Optionally store results in table:
--out assembly_overlap.tsv

"""


################################################################################

def count_bed_rows(bed_file):
    count = 0
    opener = open
    if str(bed_file).endswith(".gz"):
        import gzip
        opener = gzip.open

    with opener(bed_file, "rt") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                count += 1
    return count


################################################################################

def normalize_input_bed(in_bed, out_bed, bed_sc_thr=None, col_check=True):

    total_regions = 0
    filtered_regions = 0

    with open(in_bed) as inp, open(out_bed, "w") as out:
        for line in inp:

            if not line.strip() or line.startswith("#"):
                continue

            total_regions += 1

            fields = line.rstrip("\n").split("\t")

            if col_check and len(fields) < 6:
                assert False, f"Input BED file must have at least 6 columns (chr, start, end, name, score, strand). For BED3 support, set --unstranded and do not set --bed-sc-thr filtering (since score column is required for filtering)."

            # Optional score filter.
            if bed_sc_thr is not None:
                try:
                    score = float(fields[4])
                except ValueError:
                    assert False, f"Invalid score in line: {line.strip()}"

                if score < bed_sc_thr:
                    continue

            chr_id = benchlib.check_convert_chr_id(
                fields[0],
                id_style=1  # convert to chr1,chr2,...,chrM style.
            )

            if not chr_id:
                continue

            fields[0] = chr_id

            out.write("\t".join(fields) + "\n")
            filtered_regions += 1

    return total_regions, filtered_regions


################################################################################

def read_metadata(gene_regions_path, mode="all"):

    metadata_file = os.path.join(gene_regions_path, 'metadata.tsv')
    assert os.path.exists(metadata_file), f"Missing metadata.tsv file: {metadata_file}"

    rows = []

    with open(metadata_file, "r") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if mode != "all" and row["species"] != mode:
                continue

            # bed_path = DATA_DIR / row["file"]
            bed_path = os.path.join(gene_regions_path, row["file"])
            
            assert os.path.exists(bed_path), f"Missing gene region BED file: {bed_path}"
            
            row["bed_path"] = bed_path
            rows.append(row)

    return rows


################################################################################

def bedtools_overlap_count(query_bed, gene_bed,
                           unstranded=False):
    """
    Count input regions that overlap at least one gene region strand-specifically.

    """
    cmd = [
        "bedtools", "intersect",
        "-u",
        "-a", str(query_bed),
        "-b", str(gene_bed),
    ]
    if not unstranded:
        cmd.append("-s")

    result = subprocess.run(
        cmd,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    if not result.stdout.strip():
        return 0

    return len(result.stdout.rstrip("\n").split("\n"))


################################################################################

def main():
    # SOME GOOD ARGUMENTS (!)
    parser = argparse.ArgumentParser(
        description="Guess genome assembly of input regions BED file by strand-specific overlap with gene regions."
    )
    parser.add_argument("--in", dest="in_bed", required=True, help="Input regions BED file (>= 6 columns) to guess assembly for")
    parser.add_argument(
        "--mode",
        choices=["all", "human", "mouse"],
        default="all",
        help="Which species assemblies to test. Available assemblies: hg18 hg19 hg38 (human) and mm9 mm10 mm39 (mouse) (default: all)",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Output TSV file. By default print to STDOUT",
    )
    parser.add_argument(
        "--bed-sc-thr",
        type=float,
        default=None,
        metavar="FLOAT",
        help="Keep only BED regions with column 5 score >= FLOAT (default: no filtering)"
    )
    parser.add_argument(
        "--unstranded",
        dest="unstranded",
        default = False,
        action = "store_true",
        help="Disable strand specific overlap calculations. By default, intersecting input and gene regions have to be on same strands to count as overlapping"
    )

    args = parser.parse_args()

    in_bed = Path(args.in_bed)

    # Library path.
    benchlib_path = os.path.dirname(benchlib.__file__)
    gene_regions_path = benchlib_path + "/content/gene_regions"

    # File with gene region information from different assemblies.
    metadata_file = os.path.join(gene_regions_path, 'metadata.tsv')
    
    # Get metadata for assemblies.
    metadata_rows = read_metadata(gene_regions_path, mode=args.mode)
    
    
    """
    Filter input BED:
    1) Convert chromosome IDs to format: chr1, chr2 ...
    2) Keep only reference chromosomes.
    3) Optionally filter by BED column 5 score (--bed-sc-thr). 
    
    """

    col_check = True
    if args.bed_sc_thr is None and args.unstranded:
        # If no score filter and unstranded, we can skip column check to allow more flexible input formats (e.g. BED3).
        col_check = False

    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".filtered.chr.bed",
        delete=False,
    ) as tmp:
        tmp_bed = Path(tmp.name)

    try:
        total_regions, filtered_regions = normalize_input_bed(
            in_bed,
            tmp_bed,
            bed_sc_thr=args.bed_sc_thr,
            col_check=col_check,
        )
        
        output_rows = []

        for row in metadata_rows:
            overlapping_regions = bedtools_overlap_count(tmp_bed, row["bed_path"],
                                                         unstranded=args.unstranded)

            overlap_percent = (
                overlapping_regions / filtered_regions * 100
                if filtered_regions > 0
                else 0.0
            )

            output_rows.append({
                "assembly": row["assembly"],
                "species": row["species"],
                "source": row["source"],
                "version": row["version"],
                "total_regions": total_regions,
                "filtered_regions": filtered_regions,
                "overlapping_regions": overlapping_regions,
                "overlap_percent": f"{overlap_percent:.2f}",
                "comment": row["comment"],
                # "gene_region_file": row["file"],
            })

        # Sort assembly overlaps by overlap_percent.
        output_rows.sort(
            key=lambda x: float(x["overlap_percent"]),
            reverse=True,
        )

        fieldnames = [
            "assembly",
            "species",
            "source",
            "version",
            "total_regions",
            "filtered_regions",
            "overlapping_regions",
            "overlap_percent",
            "comment",
            # "gene_region_file",
        ]

        if args.out:
            out_handle = open(args.out, "w")
        else:
            out_handle = None

        target = out_handle if out_handle else None

        writer = csv.DictWriter(
            target or __import__("sys").stdout,
            fieldnames=fieldnames,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(output_rows)

        if out_handle:
            out_handle.close()

    finally:
        tmp_bed.unlink(missing_ok=True)


################################################################################

if __name__ == "__main__":
    main()
