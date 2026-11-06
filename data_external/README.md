# External datasets used for lifetime estimates

| file | source | content |
|---|---|---|
| `elife-34202-supp1.xlsx` | Dörrbaum, Kochen, Langer & Schuman, *eLife* 2018;7:e34202, Supplementary file 1 | SILAC protein half-lives, rat cortex/hippocampus cultures (P1-derived, 18–19 DIV), three sheets: mixed, neuron-enriched, glia-enriched. ~3,500 proteins each. Columns: index, UniProt, gene, protein, k (per day), SE, R², half-life (days). |
| `41593_2026_2420_MOESM3_ESM.xlsx` | Qiu, Zhang, Xia, Gao et al., *Nat Neurosci* 29(11):2685–2698 (2026), PMID 42675226, Supplementary Tables | Table 1: ECS-regulated ARG groups. Table 2: membership of 12,483 genes in 100 RNA kinetics modules. **Module membership only — no half-life values.** |
| `41593_2026_2420_MOESM12_ESM.xlsx` | same paper, Source Data Fig. 4 | Sheets `4b`, `4c_top_P7`, `4c_bottom_P56`, `4d`. Spatial new-to-total RNA ratios by cluster for a handful of marker genes. Notable: the `P56` sheet shows the spatial NT-seq dataset includes **adult** animals, not only P7. |

Consumed by `scripts/protein_lifetime.py` and `scripts/halflife_modules.py`.

Still missing: the per-gene RNA half-life matrix (13,591 transcripts × cell type) behind
Supplementary Fig. 19a of the Qiu paper. Not in the supplementary tables, not in the supplementary
PDF (figure images plus legends only), and not in GEO, which carries only `GSE249398_RAW.tar`
(2.0 Gb, H5AD), `GSE249403_RAW.tar` (4.6 Gb) under SuperSeries `GSE249405`. The authors' analysis
code is at `github.com/hongjie7/RNAkinetoScope` (code only, no data); `group_kinetics.py` there
documents their estimator, so a recomputation from the H5AD objects would use their own method.
