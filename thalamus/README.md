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
  cerebellum 0.71); 12% transmembrane KChIP (KChIP4a/4e) and 18% myristoylated KChIP1a/b, both far above cortex
  and hippocampus; DPP6K 13%, DPP6a 9%. No V404M-induced change survives FDR. Which cell type carries the
  TM-KChIP / Kcnd3-rich profile is what tier 1 should tell us.
