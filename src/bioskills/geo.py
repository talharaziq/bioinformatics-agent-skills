"""Optional GEO Series retrieval. Install with `pip install .[geo]`."""
from __future__ import annotations
import json
from pathlib import Path

def fetch_geo_series(accession:str, outdir:str|Path, include_matrix:bool=False)->dict:
    try: import GEOparse
    except ImportError as exc: raise RuntimeError("GEO retrieval requires optional dependency GEOparse. Install with: pip install '.[geo]'") from exc
    accession=accession.strip().upper()
    if not accession.startswith("GSE") or not accession[3:].isdigit(): raise ValueError("Expected a GEO Series accession such as GSE12345")
    out=Path(outdir); out.mkdir(parents=True,exist_ok=True)
    gse=GEOparse.get_GEO(geo=accession,destdir=str(out/"download"),silent=False)
    metadata={"accession":accession,"title":getattr(gse,"metadata",{}).get("title",[]),"summary":getattr(gse,"metadata",{}).get("summary",[]),"overall_design":getattr(gse,"metadata",{}).get("overall_design",[]),"pubmed_id":getattr(gse,"metadata",{}).get("pubmed_id",[]),"sample_count":len(getattr(gse,"gsms",{})),"platform_count":len(getattr(gse,"gpls",{}))}
    (out/"series_metadata.json").write_text(json.dumps(metadata,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    sample_path=out/"samples.tsv"
    import csv
    with sample_path.open("w",newline="",encoding="utf-8") as f:
        fields=["sample_id","title","source_name","platform_id","characteristics_json","metadata_json"]
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t"); w.writeheader()
        for gsm_id,gsm in sorted(gse.gsms.items()):
            md=getattr(gsm,"metadata",{})
            chars={k:v for k,v in md.items() if k.startswith("characteristics_ch")}
            w.writerow({"sample_id":gsm_id,"title":" | ".join(md.get("title",[])),"source_name":" | ".join(md.get("source_name_ch1",[])),"platform_id":" | ".join(md.get("platform_id",[])),"characteristics_json":json.dumps(chars,ensure_ascii=False),"metadata_json":json.dumps(md,ensure_ascii=False)})
    matrix_written=False
    if include_matrix:
        if len(getattr(gse,"gpls",{}))!=1: raise ValueError("Refusing to build one expression matrix because this series has multiple platforms. Inspect samples.tsv and process platforms separately.")
        frame=gse.pivot_samples("VALUE")
        if frame is None or frame.empty: raise ValueError("GEO series had no pivotable VALUE matrix; check the source platform and supplementary files.")
        frame.to_csv(out/"series_matrix.tsv",sep="\t",index=True)
        matrix_written=True
    return {**metadata,"samples_tsv":str(sample_path),"matrix_written":matrix_written,"output_directory":str(out.resolve())}
