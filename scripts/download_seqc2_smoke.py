#!/usr/bin/env python3
"""Download and subsample public SEQC2 ctDNA reads for a pipeline smoke test.

The source study is NCBI BioProject PRJNA677999.  Reads are streamed from ENA,
so the multi-gigabyte source FASTQs are not stored locally.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import os
import re
import sys
import tempfile
import urllib.parse
import urllib.request
from pathlib import Path

BIOPROJECT = "PRJNA677999"
ENA_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"
# Two small, paired-end AVENIO runs representing distinct SEQC2 reference
# materials. Explicit defaults keep smoke tests reproducible as ENA grows.
DEFAULT_RUNS = ("SRR13209661", "SRR13209657")  # Sample Bf and Sample Ef


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/seqc2-smoke"),
        help="Output directory (default: data/seqc2-smoke)",
    )
    parser.add_argument(
        "--reads",
        type=int,
        default=10_000,
        help="Read pairs retained per run (default: 10000)",
    )
    parser.add_argument(
        "--runs",
        nargs="+",
        default=list(DEFAULT_RUNS),
        metavar="SRR_ACCESSION",
        help="Paired-end runs to use (default: %(default)s)",
    )
    parser.add_argument("--force", action="store_true", help="Replace existing output FASTQs")
    return parser.parse_args()


def ena_metadata(run: str) -> dict[str, str]:
    query = urllib.parse.urlencode(
        {
            "accession": run,
            "result": "read_run",
            "fields": (
                "run_accession,sample_accession,sample_title,experiment_title,"
                "library_layout,fastq_ftp,fastq_md5"
            ),
            "format": "tsv",
        }
    )
    request = urllib.request.Request(
        f"{ENA_API}?{query}", headers={"User-Agent": "cfdna-demo-smoke/1.0"}
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        rows = list(csv.DictReader(io.TextIOWrapper(response, encoding="utf-8"), delimiter="\t"))
    if len(rows) != 1:
        raise RuntimeError(f"ENA returned {len(rows)} records for {run}; expected one")
    return rows[0]


def paired_urls(row: dict[str, str]) -> tuple[str, str]:
    if row["library_layout"] != "PAIRED":
        raise RuntimeError(f"{row['run_accession']} is not paired-end")
    urls = [f"https://{value}" for value in row["fastq_ftp"].split(";") if value]
    read1 = next((url for url in urls if re.search(r"_1\.f(?:ast)?q\.gz$", url)), None)
    read2 = next((url for url in urls if re.search(r"_2\.f(?:ast)?q\.gz$", url)), None)
    if not read1 or not read2:
        raise RuntimeError(f"ENA does not list paired FASTQs for {row['run_accession']}")
    return read1, read2


def normalized_read_name(header: bytes) -> bytes:
    name = header.split(maxsplit=1)[0]
    return re.sub(rb"/[12]$", b"", name)


def stream_fastq(url: str, destination: Path, count: int) -> list[bytes]:
    """Write the first count records and return their normalized read names."""
    request = urllib.request.Request(url, headers={"User-Agent": "cfdna-demo-smoke/1.0"})
    tmp_fd, tmp_name = tempfile.mkstemp(prefix=f".{destination.name}.", dir=destination.parent)
    names: list[bytes] = []
    try:
        with os.fdopen(tmp_fd, "wb") as raw_out:
            with gzip.GzipFile(fileobj=raw_out, mode="wb", mtime=0) as output:
                with urllib.request.urlopen(request, timeout=120) as response:
                    with gzip.GzipFile(fileobj=response, mode="rb") as source:
                        for record_number in range(count):
                            record = [source.readline() for _ in range(4)]
                            if not record[0]:
                                break
                            if (
                                any(not line for line in record)
                                or not record[0].startswith(b"@")
                                or not record[2].startswith(b"+")
                            ):
                                raise RuntimeError(
                                    f"Malformed FASTQ record {record_number + 1} from {url}"
                                )
                            output.writelines(record)
                            names.append(normalized_read_name(record[0][1:].rstrip()))
        os.replace(tmp_name, destination)
    except BaseException:
        try:
            os.unlink(tmp_name)
        except FileNotFoundError:
            pass
        raise
    return names


def main() -> int:
    args = parse_args()
    if args.reads < 1:
        raise SystemExit("--reads must be at least 1")
    invalid = [run for run in args.runs if not re.fullmatch(r"[SED]RR\d+", run)]
    if invalid:
        raise SystemExit(f"Invalid run accession(s): {', '.join(invalid)}")

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    samplesheet_rows: list[tuple[str, str, str]] = []
    manifest: list[dict[str, str | int]] = []

    for run in args.runs:
        row = ena_metadata(run)
        url1, url2 = paired_urls(row)
        output1 = output_dir / f"{run}_1.fastq.gz"
        output2 = output_dir / f"{run}_2.fastq.gz"
        if (output1.exists() or output2.exists()) and not args.force:
            raise RuntimeError(f"Output for {run} already exists; use --force to replace it")

        print(f"{run}: streaming {args.reads:,} read pairs from ENA", file=sys.stderr)
        names1 = stream_fastq(url1, output1, args.reads)
        names2 = stream_fastq(url2, output2, args.reads)
        if not names1 or len(names1) != len(names2) or names1 != names2:
            output1.unlink(missing_ok=True)
            output2.unlink(missing_ok=True)
            raise RuntimeError(f"Paired FASTQs for {run} contain unequal or mismatched reads")

        samplesheet_rows.append((run, str(output1), str(output2)))
        manifest.append(
            {
                "bioproject": BIOPROJECT,
                "run_accession": run,
                "sample_accession": row["sample_accession"],
                "sample_title": row["sample_title"],
                "experiment_title": row["experiment_title"],
                "read_pairs": len(names1),
                "source_fastq_1": url1,
                "source_fastq_2": url2,
            }
        )

    samplesheet = output_dir / "samplesheet.csv"
    with samplesheet.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("sample", "fastq_1", "fastq_2"))
        writer.writerows(samplesheet_rows)
    with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    print(f"Wrote {samplesheet}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
