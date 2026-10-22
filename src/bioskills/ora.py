"""Over-representation analysis (ORA) for a gene list and a local GMT collection."""
from __future__ import annotations
import csv, math
from pathlib import Path

def read_ids(path: str | Path) -> set[str]:
    p=Path(path)
    if not p.is_file(): raise FileNotFoundError(p)
    ids=set()
    for line in p.read_text(encoding="utf-8-sig").splitlines():
        value=line.strip()
        if value and not value.startswith("#"):
            ids.add(value.split("\t",1)[0].split(",",1)[0].strip())
    return ids

def read_gmt(path: str | Path) -> dict[str,set[str]]:
    terms={}
    with Path(path).open(encoding="utf-8-sig") as f:
        for no,line in enumerate(f,1):
            row=line.rstrip("\r\n").split("\t")
            if len(row)<3: raise ValueError(f"GMT line {no} must have term, description, and at least one gene")
            name=row[0].strip(); genes={x.strip() for x in row[2:] if x.strip()}
            if not name: raise ValueError(f"GMT line {no} has an empty term name")
            terms[name]=genes
    if not terms: raise ValueError("No pathways were found in the GMT file")
    return terms

def _log_choose(n:int,k:int)->float:
    if k<0 or k>n: return float("-inf")
    return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)

def hypergeom_sf(k:int, population:int, successes:int, draws:int)->float:
    """Return P[X >= k] for X~Hypergeometric(population, successes, draws)."""
    if population<=0 or successes<0 or successes>population or draws<0 or draws>population: raise ValueError("Invalid hypergeometric parameters")
    lo=max(k,0,draws-(population-successes)); hi=min(successes,draws)
    if lo>hi: return 0.0
    logs=[]
    den=_log_choose(population,draws)
    for x in range(lo,hi+1): logs.append(_log_choose(successes,x)+_log_choose(population-successes,draws-x)-den)
    mx=max(logs); result=math.exp(mx)*sum(math.exp(v-mx) for v in logs)
    return min(1.0,result)

def bh_adjust(pvalues:list[float])->list[float]:
    n=len(pvalues); order=sorted(range(n),key=lambda i:pvalues[i]); out=[1.0]*n; running=1.0
    for rank_idx in range(n-1,-1,-1):
        i=order[rank_idx]; rank=rank_idx+1
        running=min(running,pvalues[i]*n/rank); out[i]=min(1.0,running)
    return out

def run_ora(gene_list: str | Path, gmt_path: str | Path, output: str | Path, universe_path: str | Path | None = None, min_size: int = 5, max_size: int = 500) -> list[dict]:
    if min_size<1 or max_size<min_size: raise ValueError("Require 1 <= min-size <= max-size")
    query=read_ids(gene_list); terms=read_gmt(gmt_path)
    if universe_path:
        universe=read_ids(universe_path)
        outside=query-universe
        if outside: raise ValueError(f"{len(outside)} input gene(s) are not present in the stated universe; examples: {sorted(outside)[:5]}")
    else:
        universe=set(query)
        for genes in terms.values(): universe.update(genes)
    if not query: raise ValueError("Input gene list is empty")
    query &= universe; N=len(universe); n=len(query)
    if N==0 or n==0: raise ValueError("No input genes overlap with the universe")
    res=[]
    for term, all_genes in terms.items():
        members=all_genes & universe
        K=len(members)
        if K<min_size or K>max_size: continue
        overlap=query & members; k=len(overlap)
        # Keep zero-overlap terms in the multiple-testing family (p=1), even though
        # they are omitted from the presentation table below. Otherwise FDR is optimistic.
        res.append({"term":term,"overlap_count":k,"pathway_size":K,"input_size":n,"universe_size":N,"p_value":hypergeom_sf(k,N,K,n),"overlap_genes":";".join(sorted(overlap))})
    adjusted=bh_adjust([r["p_value"] for r in res])
    for r,q in zip(res,adjusted): r["fdr_bh"]=q
    res=[r for r in res if r["overlap_count"] > 0]
    res.sort(key=lambda r:(r["fdr_bh"],r["p_value"],r["term"]))
    out=Path(output); out.parent.mkdir(parents=True,exist_ok=True)
    fields=["term","overlap_count","pathway_size","input_size","universe_size","p_value","fdr_bh","overlap_genes"]
    with out.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t"); w.writeheader()
        for r in res:
            row=dict(r); row["p_value"]=f"{r['p_value']:.8g}"; row["fdr_bh"]=f"{r['fdr_bh']:.8g}"; w.writerow(row)
    return res
