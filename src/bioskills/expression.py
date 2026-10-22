"""Standard-library QC for genes-by-samples tab-separated expression matrices."""
from __future__ import annotations
import csv, json, math
from array import array
from pathlib import Path

def _pearson(x, y):
    n=len(x)
    if n < 2: return None
    mx=sum(x)/n; my=sum(y)/n
    ssx=sum((a-mx)**2 for a in x); ssy=sum((b-my)**2 for b in y)
    if ssx <= 0 or ssy <= 0: return None
    return sum((a-mx)*(b-my) for a,b in zip(x,y))/math.sqrt(ssx*ssy)

def expression_qc(matrix_path: str | Path, outdir: str | Path, metadata_path: str | Path | None = None) -> dict:
    path=Path(matrix_path); out=Path(outdir)
    if not path.is_file(): raise FileNotFoundError(f"Expression matrix not found: {path}")
    out.mkdir(parents=True, exist_ok=True)
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader=csv.reader(f, delimiter="\t")
        try: header=next(reader)
        except StopIteration: raise ValueError("Expression matrix is empty")
        if len(header)<3: raise ValueError("Expected first column gene_id plus at least two sample columns")
        sample_names=[h.strip() for h in header[1:]]
        if not header[0].strip(): raise ValueError("First column must have a header such as gene_id")
        if any(not s for s in sample_names): raise ValueError("Sample names must be non-empty")
        if len(set(sample_names)) != len(sample_names): raise ValueError("Duplicate sample names in matrix header")
        sample_values=[array('d') for _ in sample_names]
        libs=[0.0]*len(sample_names); nonzero=[0]*len(sample_names); zeros=[0]*len(sample_names)
        genes=set(); duplicates=0; fractional=0; rows=0; minimum=None; maximum=None
        for line_no,row in enumerate(reader,2):
            if not row or all(not c.strip() for c in row): continue
            if len(row)!=len(header): raise ValueError(f"Line {line_no}: expected {len(header)} columns, found {len(row)}")
            gene=row[0].strip()
            if not gene: raise ValueError(f"Line {line_no}: empty gene identifier")
            if gene in genes: duplicates+=1
            genes.add(gene); rows+=1
            for i, raw in enumerate(row[1:]):
                try: v=float(raw)
                except ValueError: raise ValueError(f"Line {line_no}, sample {sample_names[i]}: non-numeric value {raw!r}")
                if not math.isfinite(v) or v<0: raise ValueError(f"Line {line_no}, sample {sample_names[i]}: values must be finite and non-negative")
                if not v.is_integer(): fractional+=1
                sample_values[i].append(v); libs[i]+=v
                if v>0: nonzero[i]+=1
                else: zeros[i]+=1
                minimum=v if minimum is None else min(minimum,v); maximum=v if maximum is None else max(maximum,v)
    if rows==0: raise ValueError("Expression matrix contains no gene rows")
    metadata={}
    if metadata_path:
        mp=Path(metadata_path)
        if not mp.is_file(): raise FileNotFoundError(f"Sample metadata not found: {mp}")
        with mp.open("r",encoding="utf-8-sig",newline="") as f:
            reader=csv.DictReader(f,delimiter="\t")
            if not reader.fieldnames or "sample_id" not in reader.fieldnames: raise ValueError("Metadata TSV needs a sample_id column")
            for r in reader:
                sid=r.get("sample_id","").strip()
                if not sid: raise ValueError("Metadata contains an empty sample_id")
                if sid in metadata: raise ValueError(f"Duplicate sample_id in metadata: {sid}")
                metadata[sid]=r
        absent=sorted(set(sample_names)-set(metadata)); extra=sorted(set(metadata)-set(sample_names))
        if absent or extra:
            raise ValueError(f"Matrix/metadata sample mismatch; missing metadata for {absent}; metadata-only samples {extra}")
    sample_stats=[]
    for i,s in enumerate(sample_names):
        sample_stats.append({"sample_id":s,"library_sum":round(libs[i],6),"nonzero_features":nonzero[i],"zero_features":zeros[i],"zero_percent":round(100*zeros[i]/rows,4),"detected_feature_percent":round(100*nonzero[i]/rows,4),"condition":metadata.get(s,{}).get("condition","")})
    pairwise=[]
    for i in range(len(sample_names)):
        for j in range(i+1,len(sample_names)):
            r=_pearson(sample_values[i],sample_values[j])
            pairwise.append({"sample_1":sample_names[i],"sample_2":sample_names[j],"pearson_r":None if r is None else round(r,6)})
    report={"matrix":str(path),"genes_or_features":rows,"samples":len(sample_names),"duplicate_gene_rows":duplicates,"fractional_values":fractional,"value_min":minimum,"value_max":maximum,"sample_qc":sample_stats,"pairwise_pearson":pairwise,"metadata_supplied":bool(metadata_path),"warnings":[]}
    if duplicates: report["warnings"].append(f"Found {duplicates} duplicate gene IDs/rows; aggregate or resolve duplicates before downstream inference.")
    if fractional: report["warnings"].append("Fractional values detected. Do not pass normalized/fractional expression values to count-based DESeq2; use the appropriate model for the data type.")
    if not fractional: report["warnings"].append("Values are integer-like, but confirm the file contains raw, un-normalized counts before using a count model.")
    if len(sample_names)<3: report["warnings"].append("Fewer than three samples are present; check replication and design before making inferential claims.")
    with (out/"sample_qc.tsv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(sample_stats[0]),delimiter="\t"); w.writeheader(); w.writerows(sample_stats)
    with (out/"pairwise_similarity.tsv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=["sample_1","sample_2","pearson_r"],delimiter="\t"); w.writeheader(); w.writerows(pairwise)
    (out/"expression_qc.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    lines=["# Expression matrix QC report","",f"- **Features:** {rows}",f"- **Samples:** {len(sample_names)}",f"- **Duplicate feature rows:** {duplicates}",f"- **Fractional matrix values:** {fractional}","","## Per-sample summary","","| Sample | Sum of values | Detected features | Zero features (%) |","|---|---:|---:|---:|"]
    for s in sample_stats: lines.append(f"| {s['sample_id']} | {s['library_sum']:,.4g} | {s['nonzero_features']} | {s['zero_percent']:.2f}% |")
    lines += ["","## Warnings"] + [f"- {w}" for w in report["warnings"]] + ["","This report checks file integrity and coarse QC only. It does not normalize data, test differential expression, or establish biological significance.",""]
    (out/"report.md").write_text("\n".join(lines),encoding="utf-8")
    return report
