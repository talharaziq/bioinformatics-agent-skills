"""Optional Scanpy workflow for 10x directories/H5 and .h5ad inputs."""
from __future__ import annotations
import json
from pathlib import Path

def run_single_cell(input_path:str|Path,outdir:str|Path,min_genes:int=200,min_cells:int=3,max_pct_mito:float=20.0,n_top_genes:int=2000,random_state:int=0)->dict:
    try:
        import scanpy as sc
    except ImportError as exc: raise RuntimeError("Single-cell analysis requires Scanpy and its dependencies. Install with: pip install '.[single-cell]'") from exc
    inp=Path(input_path); out=Path(outdir); out.mkdir(parents=True,exist_ok=True)
    if not inp.exists(): raise FileNotFoundError(inp)
    if inp.is_file() and inp.suffix.lower()==".h5ad": adata=sc.read_h5ad(inp)
    elif inp.is_file() and inp.suffix.lower() in {".h5",".hdf5"}: adata=sc.read_10x_h5(inp)
    elif inp.is_dir(): adata=sc.read_10x_mtx(inp,var_names="gene_symbols",cache=False)
    else: raise ValueError("Input must be a .h5ad file, 10x .h5 file, or a 10x mtx directory")
    adata.var_names_make_unique(); adata.obs_names_make_unique()
    adata.var["mt"]=adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata,qc_vars=["mt"],percent_top=None,log1p=False,inplace=True)
    before=int(adata.n_obs)
    qc=adata.obs[[c for c in ["n_genes_by_counts","total_counts","pct_counts_mt"] if c in adata.obs]].copy()
    qc.index.name="cell_id"
    sc.pp.filter_cells(adata,min_genes=min_genes)
    sc.pp.filter_genes(adata,min_cells=min_cells)
    adata=adata[adata.obs.pct_counts_mt < max_pct_mito].copy()
    if adata.n_obs<4 or adata.n_vars<4: raise ValueError("Too few cells/genes remain after filtering; review thresholds and organism-specific mitochondrial gene naming")
    qc["retained_after_filtering"]=qc.index.isin(adata.obs_names)
    qc.to_csv(out/"cell_qc.tsv",sep="\t")
    # Keep the unnormalized counts in a layer, but avoid densifying/scaling all genes.
    adata.layers["counts"]=adata.X.copy()
    sc.pp.normalize_total(adata,target_sum=1e4)
    sc.pp.log1p(adata)
    hvg_n=min(n_top_genes,max(2,adata.n_vars-1))
    sc.pp.highly_variable_genes(adata,n_top_genes=hvg_n,flavor="seurat",subset=False)
    hvg=adata[:,adata.var["highly_variable"]].copy()
    if min(hvg.n_obs,hvg.n_vars) < 4:
        raise ValueError("Too few cells or highly variable genes for stable PCA; review filters/HVG settings")
    # Scale only HVGs to avoid a large dense cell-by-all-gene matrix.
    sc.pp.scale(hvg,max_value=10)
    n_comps=min(30,hvg.n_obs-1,hvg.n_vars-1)
    sc.tl.pca(hvg,n_comps=n_comps,svd_solver="arpack",random_state=random_state)
    n_pcs=min(30,hvg.obsm["X_pca"].shape[1])
    sc.pp.neighbors(hvg,n_neighbors=min(15,max(2,hvg.n_obs-1)),n_pcs=n_pcs,random_state=random_state)
    clustering="not_run"
    try:
        sc.tl.leiden(hvg,random_state=random_state)
        adata.obs["leiden"]=hvg.obs["leiden"].astype(str)
        clustering="leiden"
    except ImportError:
        pass
    sc.tl.umap(hvg,random_state=random_state)
    adata.obsm["X_pca"]=hvg.obsm["X_pca"]
    adata.obsm["X_umap"]=hvg.obsm["X_umap"]
    adata.obsp["connectivities"]=hvg.obsp["connectivities"]
    adata.obsp["distances"]=hvg.obsp["distances"]
    adata.uns["neighbors"]=hvg.uns["neighbors"]
    adata.write_h5ad(out/"processed.h5ad",compression="gzip")
    summary={"input":str(inp.resolve()),"cells_before_filtering":before,"cells_after_filtering":int(adata.n_obs),"genes_after_filtering":int(adata.n_vars),"min_genes":min_genes,"min_cells":min_cells,"max_pct_mito":max_pct_mito,"normalization":"total-count normalization to 10,000 then log1p","highly_variable_genes":int(adata.var.get("highly_variable",[]).sum()) if "highly_variable" in adata.var else None,"clustering":clustering,"outputs":["cell_qc.tsv","processed.h5ad"]}
    (out/"single_cell_summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    try:
        sc.pl.umap(adata,color=clustering if clustering=="leiden" else "total_counts",show=False,save=None)
        import matplotlib.pyplot as plt
        plt.savefig(out/"umap.png",dpi=180,bbox_inches="tight"); plt.close("all")
        summary["outputs"].append("umap.png")
        (out/"single_cell_summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    except Exception as exc:
        # Keep useful numerical output but record why the optional plot was skipped.
        summary["plot_warning"]=str(exc)
        (out/"single_cell_summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    return summary
