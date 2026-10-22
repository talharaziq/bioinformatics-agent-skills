from __future__ import annotations
import argparse, json, sys
from pathlib import Path
from .sequence import sequence_qc, write_json
from .expression import expression_qc
from .ora import run_ora
from .manifest import validate_manifest

def main(argv=None):
    parser=argparse.ArgumentParser(prog="bio-skills",description="Reproducible bioinformatics skill utilities")
    sub=parser.add_subparsers(dest="command",required=True)
    p=sub.add_parser("sequence-qc",help="Summarize FASTA/FASTQ records, lengths, GC, ambiguity, and FASTQ quality")
    p.add_argument("input"); p.add_argument("--format",choices=["auto","fasta","fastq"],default="auto"); p.add_argument("--json",dest="json_out"); p.add_argument("--records-tsv")
    p=sub.add_parser("expression-qc",help="QC a genes-by-samples expression TSV"); p.add_argument("--matrix",required=True); p.add_argument("--outdir",required=True); p.add_argument("--metadata")
    p=sub.add_parser("ora",help="Local GMT over-representation analysis with BH correction"); p.add_argument("--genes",required=True); p.add_argument("--gmt",required=True); p.add_argument("--output",required=True); p.add_argument("--universe"); p.add_argument("--min-size",type=int,default=5); p.add_argument("--max-size",type=int,default=500)
    p=sub.add_parser("validate-manifest",help="Validate RNA-seq sample sheet and input FASTQs"); p.add_argument("manifest"); p.add_argument("--check-fastq",action="store_true"); p.add_argument("--json",dest="json_out")
    p=sub.add_parser("geo-fetch",help="Download GEO Series metadata; optionally export a platform expression matrix"); p.add_argument("accession"); p.add_argument("--outdir",required=True); p.add_argument("--include-matrix",action="store_true")
    p=sub.add_parser("single-cell",help="Run optional Scanpy QC/normalization/PCA/neighbors/UMAP"); p.add_argument("input"); p.add_argument("--outdir",required=True); p.add_argument("--min-genes",type=int,default=200); p.add_argument("--min-cells",type=int,default=3); p.add_argument("--max-pct-mito",type=float,default=20.0); p.add_argument("--n-top-genes",type=int,default=2000)
    p=sub.add_parser("run-salmon",help="Run fastp, FastQC, Salmon using a sample manifest and an existing Salmon index"); p.add_argument("--manifest",required=True); p.add_argument("--index",required=True); p.add_argument("--outdir",required=True); p.add_argument("--threads",type=int,default=4); p.add_argument("--skip-fastp",action="store_true"); p.add_argument("--skip-fastqc",action="store_true")
    args=parser.parse_args(argv)
    try:
        if args.command=="sequence-qc": write_json(sequence_qc(args.input,args.format,args.records_tsv),args.json_out)
        elif args.command=="expression-qc":
            r=expression_qc(args.matrix,args.outdir,args.metadata); print(json.dumps({"outdir":str(Path(args.outdir).resolve()),"features":r["genes_or_features"],"samples":r["samples"],"warnings":r["warnings"]},indent=2))
        elif args.command=="ora":
            r=run_ora(args.genes,args.gmt,args.output,args.universe,args.min_size,args.max_size); print(json.dumps({"terms_reported":len(r),"output":str(Path(args.output).resolve()),"top_terms":[{k:v for k,v in x.items() if k in {"term","overlap_count","p_value","fdr_bh"}} for x in r[:10]]},indent=2))
        elif args.command=="validate-manifest":
            r=validate_manifest(args.manifest,args.check_fastq); rendered=json.dumps(r,indent=2)+"\n"
            if args.json_out: Path(args.json_out).write_text(rendered,encoding="utf-8")
            print(rendered,end="")
        elif args.command=="geo-fetch":
            from .geo import fetch_geo_series
            print(json.dumps(fetch_geo_series(args.accession,args.outdir,args.include_matrix),indent=2))
        elif args.command=="single-cell":
            from .singlecell import run_single_cell
            print(json.dumps(run_single_cell(args.input,args.outdir,args.min_genes,args.min_cells,args.max_pct_mito,args.n_top_genes),indent=2))
        elif args.command=="run-salmon":
            from .bulk import run_salmon
            print(json.dumps(run_salmon(args.manifest,args.index,args.outdir,args.threads,args.skip_fastp,args.skip_fastqc),indent=2))
        return 0
    except (ValueError,FileNotFoundError,RuntimeError,OSError) as exc:
        print(f"ERROR: {exc}",file=sys.stderr); return 2

if __name__=="__main__": raise SystemExit(main())
