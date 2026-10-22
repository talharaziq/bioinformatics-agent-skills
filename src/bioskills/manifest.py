"""Validate a sample sheet for paired-end/single-end FASTQ workflows."""
from __future__ import annotations
import csv, gzip
from pathlib import Path

REQUIRED={"sample_id","condition","fastq_r1","fastq_r2"}

def first_fastq_record(path:Path):
    op=gzip.open if path.name.lower().endswith((".gz",".bgz",".bgzf")) else open
    with op(path,"rt",encoding="utf-8",errors="strict") as f:
        header=f.readline().rstrip("\r\n")
        if not header.startswith("@"): raise ValueError(f"{path}: first line is not a FASTQ '@' header")
        name=header[1:].strip()
        if not name: raise ValueError(f"{path}: first FASTQ identifier is empty")
        seq=[]
        for line in f:
            line=line.rstrip("\r\n")
            if line.startswith("+"): break
            seq.append(line.strip())
        else: raise ValueError(f"{path}: first record has no '+' separator")
        seqstr="".join(seq)
        qchunks=[]; qlen=0
        for line in f:
            quality_line=line.rstrip("\r\n")
            qchunks.append(quality_line); qlen+=len(quality_line)
            if qlen>=len(seqstr): break
        if not seqstr or qlen!=len(seqstr): raise ValueError(f"{path}: first record sequence/quality length mismatch")
        if any(ord(ch)<33 or ord(ch)>126 for ch in "".join(qchunks)):
            raise ValueError(f"{path}: invalid Phred+33 quality character in first record")
        return name

def _pair_key(read_id:str)->str:
    key=read_id.split()[0] if read_id.split() else read_id
    if key.endswith("/1") or key.endswith("/2"):
        key=key[:-2]
    return key

def validate_manifest(manifest_path:str|Path, check_fastq:bool=False)->dict:
    path=Path(manifest_path)
    if not path.is_file(): raise FileNotFoundError(path)
    with path.open("r",encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f,delimiter="\t")
        header=reader.fieldnames or []
        if len(header)!=len(set(header)): raise ValueError("Manifest has duplicate column names")
        fields=set(header)
        missing=REQUIRED-fields
        if missing: raise ValueError(f"Manifest missing columns: {', '.join(sorted(missing))}. Required: {', '.join(sorted(REQUIRED))}")
        metadata_columns=[c for c in header if c not in (REQUIRED | {"layout"})]
        rows=list(reader)
    if not rows: raise ValueError("Manifest has no samples")
    ids=set(); conditions=set(); details=[]; base=path.resolve().parent
    for line_no,row in enumerate(rows,2):
        sid=(row.get("sample_id") or "").strip(); cond=(row.get("condition") or "").strip()
        r1=(row.get("fastq_r1") or "").strip(); r2=(row.get("fastq_r2") or "").strip()
        if not sid: raise ValueError(f"Line {line_no}: sample_id is empty")
        if any(c in sid for c in "/\\\t\n\r") or sid in {".",".."}: raise ValueError(f"Line {line_no}: sample_id is unsafe for directory names: {sid!r}")
        if sid in ids: raise ValueError(f"Line {line_no}: duplicate sample_id {sid!r}")
        ids.add(sid)
        if not cond: raise ValueError(f"Line {line_no}: condition is empty for {sid}")
        conditions.add(cond)
        if not r1: raise ValueError(f"Line {line_no}: fastq_r1 is required for {sid}")
        p1=Path(r1); p1=(base/p1).resolve() if not p1.is_absolute() else p1.resolve()
        p2=None
        if r2:
            p2=Path(r2); p2=(base/p2).resolve() if not p2.is_absolute() else p2.resolve()
            if p1==p2: raise ValueError(f"Sample {sid}: R1 and R2 point to the same file")
        first_ids={}
        for label,p in [("R1",p1)]+([("R2",p2)] if p2 else []):
            if not p.is_file(): raise FileNotFoundError(f"Sample {sid} {label} file not found: {p}")
            if p.stat().st_size==0: raise ValueError(f"Sample {sid} {label} file is empty: {p}")
            if check_fastq: first_ids[label]=first_fastq_record(p)
        if check_fastq and p2 and _pair_key(first_ids["R1"]) != _pair_key(first_ids["R2"]):
            raise ValueError(f"Sample {sid}: first read identifiers differ between R1 ({first_ids['R1']!r}) and R2 ({first_ids['R2']!r})")
        extra={column:(row.get(column) or "").strip() for column in metadata_columns}
        details.append({"sample_id":sid,"condition":cond,"layout":"paired-end" if p2 else "single-end","fastq_r1":str(p1),"fastq_r2":str(p2) if p2 else "",**extra})
    if len(conditions)<2: raise ValueError("At least two condition groups are required for a comparison")
    counts={c:sum(r["condition"]==c for r in details) for c in sorted(conditions)}
    warnings=[]
    for c,n in counts.items():
        if n<2: warnings.append(f"Condition {c!r} has {n} sample(s); biological replication is limited and inference may be unreliable.")
    return {"manifest":str(path.resolve()),"n_samples":len(details),"conditions":counts,"metadata_columns":metadata_columns,"samples":details,"warnings":warnings}
