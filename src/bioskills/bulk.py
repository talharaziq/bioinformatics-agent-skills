"""Small, auditable orchestrator for fastp + FastQC + Salmon."""
from __future__ import annotations
import csv, shutil, subprocess
from pathlib import Path
from .manifest import validate_manifest

def run_salmon(manifest_path:str|Path,index_path:str|Path,outdir:str|Path,threads:int=4,skip_fastp:bool=False,skip_fastqc:bool=False)->dict:
    if threads<1: raise ValueError("threads must be >= 1")
    index=Path(index_path).resolve(); out=Path(outdir).resolve(); out.mkdir(parents=True,exist_ok=True)
    if not index.exists(): raise FileNotFoundError(f"Salmon index not found: {index}")
    for tool in ["salmon"]+([] if skip_fastqc else ["fastqc"])+([] if skip_fastp else ["fastp"]):
        if not shutil.which(tool): raise RuntimeError(f"Required executable '{tool}' was not found on PATH. Install it or choose the corresponding skip option.")
    m=validate_manifest(manifest_path)
    fastqc_dir=out/"fastqc"; salmon_dir=out/"salmon"; trim_dir=out/"trimmed"; logs=out/"logs"
    for p in [fastqc_dir,salmon_dir,trim_dir,logs]: p.mkdir(parents=True,exist_ok=True)
    existing=[s["sample_id"] for s in m["samples"] if (salmon_dir/s["sample_id"]).exists()]
    if existing:
        raise FileExistsError("Salmon output already exists for sample(s): " + ", ".join(existing) + ". Move/delete those folders or choose a new output directory.")
    for s in m["samples"]:
        sid=s["sample_id"]; r1=Path(s["fastq_r1"]); r2=Path(s["fastq_r2"]) if s["fastq_r2"] else None
        quant=salmon_dir/sid
        if not skip_fastqc:
            per_sample_fastqc=fastqc_dir/sid; per_sample_fastqc.mkdir(parents=True,exist_ok=True)
            cmd=["fastqc","--threads",str(threads),"--outdir",str(per_sample_fastqc),str(r1)]+([str(r2)] if r2 else [])
            _run(cmd,logs/f"{sid}.fastqc.log")
        if skip_fastp:
            q1,q2=r1,r2
        else:
            q1=trim_dir/f"{sid}_R1.fastq.gz"; q2=trim_dir/f"{sid}_R2.fastq.gz" if r2 else None
            cmd=["fastp","--in1",str(r1),"--out1",str(q1),"--thread",str(threads),"--html",str(trim_dir/f"{sid}.fastp.html"),"--json",str(trim_dir/f"{sid}.fastp.json")]
            if r2: cmd += ["--in2",str(r2),"--out2",str(q2)]
            _run(cmd,logs/f"{sid}.fastp.log")
        cmd=["salmon","quant","-i",str(index),"-l","A","-p",str(threads),"--validateMappings","-o",str(quant)]
        cmd += (["-1",str(q1),"-2",str(q2)] if q2 else ["-r",str(q1)])
        _run(cmd,logs/f"{sid}.salmon.log")
    output_fields=["sample_id","condition","layout","fastq_r1","fastq_r2"]+m["metadata_columns"]
    with (out/"run_manifest.tsv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=output_fields,delimiter="\t"); w.writeheader(); w.writerows(m["samples"])
    return {"samples":m["n_samples"],"conditions":m["conditions"],"output_directory":str(out),"next_step":"Run scripts/tximport_deseq2.R with a transcript_id-to-gene_id mapping to test differential expression."}

def _run(cmd:list[str],log:Path):
    with log.open("w",encoding="utf-8") as fh:
        fh.write("COMMAND: "+" ".join(cmd)+"\n")
        proc=subprocess.run(cmd,stdout=fh,stderr=subprocess.STDOUT,text=True)
    if proc.returncode:
        raise RuntimeError(f"Command failed (exit {proc.returncode}); inspect log: {log}")
