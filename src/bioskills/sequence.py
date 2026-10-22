"""Streaming FASTA/FASTQ quality summaries using only the Python standard library."""
from __future__ import annotations
import csv, gzip, json
from collections import Counter
from pathlib import Path
from typing import Iterator, TextIO

CANONICAL_NUCLEOTIDES = set("ACGTU")
IUPAC_NUCLEOTIDES = set("ACGTUNRYKMSWBDHV.-")

class SequenceFormatError(ValueError):
    """Raised when a sequence file does not conform to FASTA or FASTQ."""

def _open_text(path: Path) -> TextIO:
    if str(path).lower().endswith((".gz", ".bgz", ".bgzf")):
        return gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="")
    return path.open("rt", encoding="utf-8", errors="strict", newline="")

def _lines(handle: TextIO) -> Iterator[tuple[int, str]]:
    for i, line in enumerate(handle, 1):
        yield i, line.rstrip("\r\n")

def fasta_records(handle: TextIO) -> Iterator[tuple[str, str, None]]:
    name = None
    chunks: list[str] = []
    for lineno, line in _lines(handle):
        if not line.strip():
            continue
        if line.startswith(">"):
            if name is not None:
                seq = "".join(chunks).upper()
                if not seq:
                    raise SequenceFormatError(f"Empty FASTA sequence before line {lineno} ({name})")
                yield name, seq, None
            name = line[1:].strip()
            if not name:
                raise SequenceFormatError(f"Missing FASTA identifier on line {lineno}")
            chunks = []
        else:
            if name is None:
                raise SequenceFormatError(f"Sequence data before first FASTA header on line {lineno}")
            chunks.append("".join(line.split()))
    if name is not None:
        seq = "".join(chunks).upper()
        if not seq:
            raise SequenceFormatError(f"Empty FASTA sequence ({name})")
        yield name, seq, None
    elif not chunks:
        raise SequenceFormatError("No FASTA records found")

def fastq_records(handle: TextIO) -> Iterator[tuple[str, str, str]]:
    """Parse FASTQ records, including wrapped sequence/quality lines."""
    it = iter(_lines(handle))
    while True:
        try:
            lineno, header = next(it)
        except StopIteration:
            return
        if not header:
            continue
        if not header.startswith("@"):
            raise SequenceFormatError(f"Expected FASTQ header '@' at line {lineno}")
        name = header[1:].strip()
        if not name:
            raise SequenceFormatError(f"Missing FASTQ identifier at line {lineno}")
        seq_chunks: list[str] = []
        for seq_line_no, line in it:
            if line.startswith("+"):
                break
            if not line:
                continue
            seq_chunks.append("".join(line.split()))
        else:
            raise SequenceFormatError(f"Missing '+' separator for record {name}")
        seq = "".join(seq_chunks).upper()
        if not seq:
            raise SequenceFormatError(f"Empty FASTQ sequence for record {name}")
        qual_chunks: list[str] = []
        qual_len = 0
        for qline_no, line in it:
            qual_chunks.append(line)
            qual_len += len(line)
            if qual_len >= len(seq):
                break
        if qual_len != len(seq):
            raise SequenceFormatError(f"Quality length {qual_len} != sequence length {len(seq)} for {name}")
        quality = "".join(qual_chunks)
        if any(ord(ch) < 33 or ord(ch) > 126 for ch in quality):
            raise SequenceFormatError(f"Quality score contains a character outside printable Phred+33 range for {name}")
        yield name, seq, quality

def _median_from_hist(hist: Counter[int], n: int) -> float:
    if n == 0:
        return 0.0
    ranks = [(n + 1) // 2, (n + 2) // 2]
    found = []
    seen = 0
    for length in sorted(hist):
        seen += hist[length]
        while ranks and seen >= ranks[0]:
            found.append(length); ranks.pop(0)
        if not ranks:
            break
    return sum(found) / len(found)

def sequence_qc(input_path: str | Path, fmt: str = "auto", records_tsv: str | Path | None = None) -> dict:
    path = Path(input_path)
    if not path.is_file():
        raise FileNotFoundError(f"Input file does not exist: {path}")
    selected = fmt.lower()
    if selected == "auto":
        low = path.name.lower().removesuffix(".gz").removesuffix(".bgz")
        if low.endswith((".fa", ".fasta", ".fna", ".ffn", ".faa")):
            selected = "fasta"
        elif low.endswith((".fq", ".fastq")):
            selected = "fastq"
        else:
            with _open_text(path) as h:
                first = next((ln for ln in h if ln.strip()), "")
            selected = "fasta" if first.startswith(">") else "fastq" if first.startswith("@") else ""
    if selected not in {"fasta", "fastq"}:
        raise SequenceFormatError("Could not infer format; specify --format fasta or fastq")

    hist: Counter[int] = Counter()
    n_records = total_bases = gc = ambiguous = quality_sum = quality_bases = 0
    base_counts: Counter[str] = Counter()
    out_handle = None
    writer = None
    if records_tsv:
        outpath = Path(records_tsv); outpath.parent.mkdir(parents=True, exist_ok=True)
        out_handle = outpath.open("w", newline="", encoding="utf-8")
        writer = csv.DictWriter(out_handle, fieldnames=["id", "length", "gc_percent", "ambiguous_bases", "mean_phred33"], delimiter="\t")
        writer.writeheader()
    try:
        with _open_text(path) as handle:
            records = fasta_records(handle) if selected == "fasta" else fastq_records(handle)
            for name, seq, quality in records:
                invalid = sorted(set(seq) - IUPAC_NUCLEOTIDES)
                if invalid:
                    raise SequenceFormatError(f"Record {name!r} contains non-nucleotide symbol(s) {invalid}; this tool is intended for nucleotide sequences, not protein FASTA")
                length = len(seq); n_records += 1; total_bases += length; hist[length] += 1
                base_counts.update(seq)
                gc_i = seq.count("G") + seq.count("C"); gc += gc_i
                amb_i = sum(1 for c in seq if c not in CANONICAL_NUCLEOTIDES)
                ambiguous += amb_i
                qmean = ""
                if quality is not None:
                    qsum = sum(max(0, ord(ch) - 33) for ch in quality)
                    quality_sum += qsum; quality_bases += len(quality)
                    qmean = round(qsum / len(quality), 4) if quality else 0
                if writer:
                    writer.writerow({"id": name, "length": length, "gc_percent": round(100 * gc_i / length, 4), "ambiguous_bases": amb_i, "mean_phred33": qmean})
    finally:
        if out_handle:
            out_handle.close()
    if not n_records:
        raise SequenceFormatError("No sequence records found")
    result = {
        "input": str(path), "format": selected, "records": n_records, "total_bases": total_bases,
        "mean_length": round(total_bases / n_records, 4), "median_length": _median_from_hist(hist, n_records),
        "min_length": min(hist), "max_length": max(hist), "gc_percent": round(100 * gc / total_bases, 4),
        "ambiguous_bases": ambiguous, "ambiguous_percent": round(100 * ambiguous / total_bases, 4),
        "base_counts": {k: base_counts.get(k, 0) for k in sorted(set("ACGTUNRYKMSWBDHV.-") | set(base_counts))},
        "length_distribution": {str(k): hist[k] for k in sorted(hist)},
    }
    if selected == "fastq":
        result["mean_phred33_quality"] = round(quality_sum / quality_bases, 4) if quality_bases else None
        result["quality_encoding_assumption"] = "Phred+33 (verify instrument/pipeline provenance for older FASTQ files)"
    return result

def write_json(result: dict, output: str | Path | None) -> None:
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if output:
        p = Path(output); p.parent.mkdir(parents=True, exist_ok=True); p.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
