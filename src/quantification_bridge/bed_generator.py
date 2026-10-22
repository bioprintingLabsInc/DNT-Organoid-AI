"""Deterministic Ensembl Release 112 GTF to BED12 converter for RSeQC."""

import gzip
import hashlib
from collections import defaultdict
from pathlib import Path


def sha256_file(path: str | Path) -> str:
    """Calculate SHA-256 hash of a file."""
    sha = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            sha.update(chunk)
    return sha.hexdigest()


def generate_bed12_from_gtf(
    gtf_path: str | Path,
    output_bed_path: str | Path,
) -> dict[str, str | int]:
    """Convert Ensembl GTF to standard BED12 format deterministically for RSeQC.

    Output format:
    1. chrom: chromosome name without 'chr' prefix (matching Ensembl FASTA and STAR BAMs)
    2. chromStart: 0-based start
    3. chromEnd: 1-based end
    4. name: transcript_id
    5. score: 0
    6. strand: '+' or '-'
    7. thickStart: CDS start or chromStart
    8. thickEnd: CDS end or chromEnd
    9. itemRgb: 0
    10. blockCount: number of exons
    11. blockSizes: comma-separated list of exon lengths
    12. blockStarts: comma-separated list of exon start offsets relative to chromStart
    """
    gtf_p = Path(gtf_path)
    out_p = Path(output_bed_path)

    if not gtf_p.exists():
        raise FileNotFoundError(f"Input GTF file not found: {gtf_p}")

    gtf_sha = sha256_file(gtf_p)

    # Data structures: transcript_id -> {chrom, strand, exons: list of (start, end), cds: list of (start, end)}
    transcripts = {}
    
    open_func = gzip.open if str(gtf_p).endswith(".gz") else open
    with open_func(gtf_p, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            
            feature = parts[2]
            if feature not in ("exon", "CDS"):
                continue

            chrom = parts[0]
            start = int(parts[3]) - 1  # 1-based GTF to 0-based BED
            end = int(parts[4])
            strand = parts[6]
            attrs_str = parts[8]

            # Parse transcript_id
            t_id = None
            for item in attrs_str.strip(";").split(";"):
                item = item.strip()
                if item.startswith("transcript_id "):
                    t_id = item.split(" ", 1)[1].strip('"')
                    break

            if not t_id:
                continue

            if t_id not in transcripts:
                transcripts[t_id] = {
                    "chrom": chrom,
                    "strand": strand,
                    "exons": [],
                    "cds": [],
                }

            if feature == "exon":
                transcripts[t_id]["exons"].append((start, end))
            elif feature == "CDS":
                transcripts[t_id]["cds"].append((start, end))

    # Build BED12 records
    bed_records = []
    for t_id, data in transcripts.items():
        exons = sorted(set(data["exons"]), key=lambda x: x[0])
        if not exons:
            continue

        chrom = data["chrom"]
        strand = data["strand"]
        tx_start = exons[0][0]
        tx_end = exons[-1][1]

        # CDS thick bounds
        if data["cds"]:
            cds_sorted = sorted(set(data["cds"]), key=lambda x: x[0])
            thick_start = cds_sorted[0][0]
            thick_end = cds_sorted[-1][1]
        else:
            thick_start = tx_start
            thick_end = tx_end

        block_count = len(exons)
        block_sizes = [str(e - s) for s, e in exons]
        block_starts = [str(s - tx_start) for s, e in exons]

        record = (
            chrom,
            tx_start,
            tx_end,
            t_id,
            0,
            strand,
            thick_start,
            thick_end,
            0,
            block_count,
            ",".join(block_sizes) + ",",
            ",".join(block_starts) + ",",
        )
        bed_records.append(record)

    # Sort deterministically: by chrom (alphanumeric), tx_start, t_id
    def chrom_sort_key(c: str) -> tuple[int, str]:
        if c.isdigit():
            return (0, f"{int(c):04d}")
        return (1, c)

    bed_records.sort(key=lambda r: (chrom_sort_key(r[0]), r[1], r[2], r[3]))

    # Write output BED12
    out_p.parent.mkdir(parents=True, exist_ok=True)
    with out_p.open("w", encoding="utf-8", newline="\n") as f:
        for rec in bed_records:
            f.write(f"{rec[0]}\t{rec[1]}\t{rec[2]}\t{rec[3]}\t{rec[4]}\t{rec[5]}\t{rec[6]}\t{rec[7]}\t{rec[8]}\t{rec[9]}\t{rec[10]}\t{rec[11]}\n")

    out_sha = sha256_file(out_p)

    return {
        "status": "BED12_GENERATION_SUCCESS",
        "gtf_path": str(gtf_p),
        "gtf_sha256": gtf_sha,
        "output_bed_path": str(out_p),
        "output_bed_sha256": out_sha,
        "total_transcripts_written": len(bed_records),
    }
