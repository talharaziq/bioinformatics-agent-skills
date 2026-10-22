---
name: pathway-enrichment
description: Test whether a submitted gene list is over-represented in local GMT pathways or gene sets. Use when asked to interpret a differential-expression hit list, run GO/pathway ORA, rank enriched terms, or calculate hypergeometric p-values and Benjamini-Hochberg FDR.
---

# Pathway over-representation analysis

## Procedure

1. Identify the source contrast and select the significant-gene list using a declared threshold (for example adjusted p-value < 0.05 plus a direction/effect-size rule when relevant). Do not choose a threshold only to make a pathway appear significant.
2. Confirm identifier namespace and version for the hit list, universe, and GMT file. Exact matching is used; the script does not map aliases, convert Ensembl IDs, or remove version suffixes automatically.
3. Supply the universe of genes that could have entered the hit list. For differential expression, this usually means genes eligible for testing after prefiltering, not every gene in the genome and not only the detected hits.
4. Run ORA using the desired GMT collection. The script calculates the hypergeometric upper-tail probability for each qualifying set and Benjamini–Hochberg adjusted FDR across the sets actually tested.
5. Inspect the overlap genes and pathway sizes for every top term; report input size, universe size, collection/version, gene identifiers, size filters, p-value, and FDR. Check redundancy between related terms.
6. Treat terms as descriptive enrichment, not evidence of causality or pathway activation. ORA is not a ranked-list test and does not preserve the sign of gene effects. Analyze up- and downregulated lists separately if that matches the question.

## Run

```bash
bio-skills ora --genes significant_genes.txt --gmt pathways.gmt --universe tested_genes.txt --output results/pathway_ora.tsv --min-size 5 --max-size 500
```

A tiny synthetic GMT, gene list, and universe are provided in `examples/pathways/`. The helper script is `skills/pathway-enrichment/scripts/run_ora.py`.

## Method details

For universe size N, pathway size K, hit-list size n, and observed overlap k, the one-sided p-value is P(X ≥ k) for a hypergeometric random variable. Gene-set sizes are computed after intersecting each set with the stated universe. Terms with no overlap are omitted from the presentation table but still count in the multiple-testing family with p = 1; qualifying terms are corrected with the Benjamini–Hochberg procedure.

## Quality gates and pitfalls

- A missing/incorrect universe changes the null model and can alter results substantially.
- Identifier mismatch can produce empty or severely reduced overlaps; investigate before reporting “no enrichment.”
- Redundant gene sets yield correlated terms; FDR correction does not make pathways independent.
- A pathway database's name, release/date, organism, and license must be recorded. The repository does not bundle a large proprietary/copyright-restricted pathway database.
- Very small gene lists and low-complexity sets have limited power; lack of significance is not proof of absence.
