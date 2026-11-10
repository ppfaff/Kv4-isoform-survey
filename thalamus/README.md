# Thalamus extension — Kv4 complex across the thalamocortical loop (2026-11-09)

Goal: which thalamic neurons (reticular nucleus, first-order and higher-order relay, intralaminar) express
Kv4.2 at all, with which Kv4 paralogs and which KChIP/DPP partners, as the cellular substrate for the
spike-wave discharges in Kcnd2-V404M/+ mice. Three tiers:

| tier | data | resolution | status |
|---|---|---|---|
| 1 | ABC Atlas 10x whole-brain (thalamus dissections) + THALMANAC MERFISH thalamus subset | gene level, cell type and CCFv3 nucleus | scripts here; run on a lab machine |
| 2 | ThalamoSeq (Phillips 2019, GSE133911 pooled SMART-seq by nucleus; GSE133912 single cells) + a full-length TRN dataset | N-terminal variant level | next |
| 3 | in-house V404M vs WT bulk (thalamus libraries) | variant level, whole tissue | done: `scripts/v404m_thalamus_view.py` |

THALMANAC (github.com/AllenNeuralDynamics/thalmanac; package `abc-merfish-analysis`) gives the TH+ZI spatial
subset of the ABC Atlas MERFISH data with CCFv3 nucleus labels, a cluster-to-nucleus annotation table
(vendored here as `annotations_c2n_combined.csv`) and the thalamic slice of the CCN20230722 taxonomy
(`taxonomy_CCN20230722_th_zi_only.csv`). The MERFISH panel is 500 cell-typing genes and does not measure the
Kv4 complex directly, so expression is read from the 10x data per cluster and painted onto the MERFISH
cells by cluster label (what Allen's 'imputed' release does, restricted to our panel and deterministic).

## Files
- `th_common.py` — gene panel (Kv4 complex, burst machinery, markers), thalamocortical-loop grouping of nuclei,
  cluster→nucleus lookup, parquet/csv helpers
- `extract_th_10x.py` — pull the panel for thalamic neurons from WMB-10Xv3-TH / WMB-10Xv2-TH (`--extra-matrices`
  also sweeps HY/MB for RT-ZI cells dissected outside TH); writes `tables/th10x_cells.parquet`
- `summarize_th_10x.py` — per cluster / nucleus / subclass / loop element tables and three figures
- `merfish_nucleus_map.py` — THALMANAC MERFISH subset painted with 10x cluster means; per-nucleus table,
  painted sections, measured-vs-painted check for any panel gene the MERFISH panel does carry

## Run (WSL on PAULP-ALIEN or PFAFF-LAMBDA; ~10 GB download for the two TH matrices, ~30 GB with MERFISH)
```bash
micromamba create -n abc python=3.11 pandas numpy scipy matplotlib anndata pyarrow nibabel networkx -c conda-forge
micromamba activate abc
pip install git+https://github.com/AllenInstitute/abc_atlas_access
pip install git+https://github.com/AllenNeuralDynamics/abc-merfish-analysis
cd ~ && git clone https://github.com/ppfaff/Kv4-isoform-survey   # or git pull
mkdir -p ~/kv4_thalamus/tables && cd ~/kv4_thalamus
CACHE=/mnt/d/abc_atlas          # any disk with ~50 GB free
python3 ~/Kv4-isoform-survey/thalamus/extract_th_10x.py --cache $CACHE --out tables            # TH dissections only
# python3 ~/Kv4-isoform-survey/thalamus/extract_th_10x.py --cache $CACHE --out tables --extra-matrices   # + RT-ZI from HY/MB
python3 ~/Kv4-isoform-survey/thalamus/summarize_th_10x.py tables .
python3 ~/Kv4-isoform-survey/thalamus/merfish_nucleus_map.py --cache $CACHE --cluster-means tables/th10x_cluster_means.csv --out .
```
Copy `tables/th10x_*.csv`, `tables/merfish_nucleus_kv4.csv` and `figures/*.png` back into the repo.
The per-cell tables (`*_cells.parquet`) stay on the machine (same convention as the cortical per-cell tables).

Smoke-tested on synthetic cells (taxonomy labels real, expression random); the 10x `log2` matrices are
log2(CPM+1), which is what the scripts assume when they back-transform to CPM for averaging.

## Notes for interpretation
- Loop groups (th_common.LOOP_GROUPS): TRN (RT); first-order relay VPM/VPL/LGd/MG; higher-order PO/LP/MD;
  motor VAL/VM; intralaminar/midline CM/CL/PCN/PF/RE…; anterior AD/AV/AM/LD/PVT; epithalamus; ZI/LGv/IGL.
- Clusters of `151 TH Prkcd Grin2c Glut` (the big sensory/motor relay subclass) are resolved to nuclei only
  through the THALMANAC annotation table; unannotated ones are reported as "relay (nucleus unresolved)".
- Tier 3 result (bulk): thalamus Kcnd2 share of Kv4 = 0.50 (lowest of the four tissues; cortex 0.64,
  cerebellum 0.71), Kcnd3 0.45; no V404M-induced gene-level change. **The bulk N-terminal variant fractions
  (panels C-E of thalamus_bulk_kv4_complex.png) are NOT reliable**: the 3'-biased polyA libraries collapse
  5'-distinguished variants into one sink per gene. Shown by simulation for DPP (lab_db write-up) and now
  empirically for KChIP against the Tasic 2016 pseudo-bulk (scripts/bulk_vs_tasic_variants.py,
  tables/bulk_vs_tasic_variant_calibration.csv): KChIP2x 0.77 -> 0.001, KChIP1c 0.19 -> 0.01, DPP6a 0.53 -> 0.01.
  Only Kcnip3 (3a / cys-rich), Kv4.3S/L and the gene-level numbers survive. Thalamic variant composition
  therefore has to come from tier 2 (full-length data), not from the bulk.

## Tier-1 results (2026-11-10), ABC Atlas 10x, 183,553 thalamic neurons — tables/th10x_*.csv, figures/th10x_kv4_calibrated.png
Cross-gene ratios are calibrated for gene length by anchoring each gene's whole-thalamus 10x mean to the bulk exon CPM
(thalamus/calibrate_th10x.py; within-gene contrasts between cell types are the raw 10x values).
- **TRN is the most Kv4.2-dependent thalamic population**: Kcnd3 ~12 CPM vs 200-340 in every relay group (raw), calibrated
  Kcnd2 share 0.77, Kcnd2 in 96% of cells; partners KChIP1 (57%) + KChIP2 (36%), DPP10-dominant, KChIP:Kv4 ~3.4,
  DPP:Kv4 ~8.9 - the cortical interneuron signature. Slc32a1 in 93% of RT cells (Viaat-Cre is a clean TRN driver);
  Camk2a transcript in 83% of RT cells at half the relay level.
- **Hartley 2024 subnetworks** (thalamus/trn_subnetworks.py, per cell; tables/th10x_trn_core_vs_shell.csv): 415 Spp1+ core
  cells (4.6%), concentrated in supertype Gnb3_5 (cluster 1604 is 79% Spp1+). Core vs rest: Kcnd3 7 vs 38 (0.5% of Kcnd2),
  Kcnip1 830 vs 410, Kcnip4 109 vs 450, Dpp10 3,020 vs 2,250, Cacna1h detected 59% vs 19%; differences hold within
  Pvalb-high cells. The FO-projecting core is the purest Kv4.2-only population in the dataset, with a Kv4.2-KChIP1-DPP10
  complex and the strongest T-type expression. Spp1-Cre x cKI is the targeted test.
- **Relay neurons co-express Kcnd2/Kcnd3** with a gradient: first-order sensory and anterior/motor most Kv4.3-balanced
  (LGd 0.31, VPM/VPL 0.30, LD 0.28, AD 0.34, AV 0.42); higher-order/midline/intralaminar Kcnd2-dominant (LP/SGN 0.70,
  PVT/IMD 0.69, Xi 0.70, CM 0.62). Partners KChIP3 (first-order) or KChIP4 (intralaminar/anterior/motor) with DPP6 60-75%.
  Kcnd1 negligible everywhere. Medial habenula: KChIP1 88%.
- Tier-2 variant targets: TRN = KChIP1a/b/c/e, KChIP2a/2x, DPP10a/d (cortical Pvalb interneuron variants as first proxy);
  relay = KChIP3a/3x, KChIP4a/b/bL/d/e, DPP6a/K/S.

## Tier-2 proxy (2026-11-10): TRN variants predicted from cortical interneurons (Tasic 2016, already analysed)
TRN shares lineage and gene-level complex with cortical Pvalb/Sst interneurons; their variant fractions are the
working guess for TRN until full-length TRN data exist (tables/ctx_variant_frac_of_gene_by_type.csv, GABA types).

| gene | Pvalb cells (proxy for Spp1+ core) | Sst cells (proxy for shell) |
| --- | --- | --- |
| Kcnip1 | 1c 32%, 1a 23%, 1b 18%, 1e 6% (TM), new-short 20% | 1b 41%, 1a 27%, 1c 12%, 1e 0% |
| Kcnip2 | 2x 78% (TM), 2a 18% | 2x 71%, 2a 25% (but 6x less Kcnip2 than Pvalb) |
| Dpp10 | 10d 86%, 10c 13%, 10a 1% | 10d 79%, 10a 13%, 10c 8% |
| Dpp6 | 6a 56%, 6S 29%, 6K 9% | 6a 54%, 6K 27%, 6S 15% |

Predicted core-TRN complex: Kv4.2 + KChIP1 (myristoyl/palmitoyl mix, some 1e TM) + KChIP2x + DPP10d; no Kv4.3, little KChIP4.
Measurement requires full-length RNA-seq of TRN (ThalamoSeq likely excludes RT; a Viaat-Cre/Spp1-Cre;Ai14 sort with
SMART-seq would give both the variants and the cKI recombination efficiency in the same cells).
