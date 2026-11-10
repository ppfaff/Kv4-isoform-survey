# The Kv4 complex across the thalamocortical loop: which thalamic neurons depend on Kv4.2, and with which auxiliary subunits

Nov 10, 2026 · Paul Pfaffinger (analysis with Claude) · PROVISIONAL — relay-neuron variant section pending ThalamoSeq run

## Summary

To ask which thalamic neurons are substrates for the spike-wave discharges of Kv4.2-V404M/+ mice, the Allen Brain Cell Atlas 10x whole-brain data were subset to 183,553 thalamic neurons (TH and HY/MB dissections; the THALMANAC taxonomy slice and cluster-to-nucleus annotations), summarized per thalamocortical-loop element and per CCFv3 nucleus, and anchored for gene-length bias to the lab's bulk thalamus exon counts.

**The reticular nucleus is the most Kv4.2-dependent neuronal population in the thalamus.** RT neurons express Kcnd3 at ~12 CPM against 200–340 in every relay group (a 10–30-fold within-gene deficit), so their A-current is built on Kv4.2 with essentially no Kv4.3 to buffer a dominant-negative mutant; calibrated Kcnd2 share is 0.77, with Kcnd2 detected in 96% of RT cells. Their auxiliary complement is the cortical-interneuron one — KChIP1 (57%) + KChIP2 (36%), DPP10-dominant, with auxiliary subunits in large excess (KChIP:Kv4 ≈ 3.4, DPP:Kv4 ≈ 8.9) — and shares nothing with relay neurons, which use KChIP3/KChIP4 + DPP6. Within RT, the Spp1⁺ "core" cells of Hartley et al. (2024; first-order-projecting, Spp1-Cre-accessible) are the extreme case: Kcnd3 at 0.5% of Kcnd2, KChIP1 twice and KChIP4 a quarter of the rest of RT, DPP10 highest, and the strongest T-type channel expression (Cacna1h detected in 59% vs 19%). They are concentrated in the Allen supertype RT-ZI Gnb3_5.

**Relay neurons co-express Kcnd2 and Kcnd3, with a gradient.** First-order sensory and anterior/motor relays are the most Kv4.3-balanced (LGd 0.31, VPM/VPL 0.30, LD 0.28, AD 0.34, AV 0.42, VAL 0.41 Kcnd2 share); higher-order, midline and intralaminar nuclei are Kcnd2-dominant (LP/SGN 0.70, PVT/IMD 0.69, Xi 0.70, CM 0.62, IAD 0.60). Partners are KChIP3 (first-order sensory) or KChIP4 (intralaminar, anterior, motor) with DPP6 at 60–75%. Kcnd1 is negligible everywhere.

**Cre coverage from the same cells**: Slc32a1 (Viaat) is detected in 93% of RT cells and 1–2% of relay cells, so Viaat-Cre is a clean TRN driver within thalamus; Camk2a transcript is present in 83% of RT cells at half the relay level, so a CaMKIIα-Cre transgene's TRN expression has to be checked with the reporter, not assumed. Vglut2-Cre covers relay neurons (Slc17a6 in 100%) and not RT (16%).

Read with the lab's conditional results — Viaat-Cre × cKI reproduces the spike-wave discharges and shows ΔFosB in TRN and upper-layer cortex but not in dentate/CA1 — the thalamic side of the loop phenotype has an obvious cellular substrate: V404M in a Kv4.2-only, KChIP1/KChIP2/DPP10-equipped TRN neuron with strong burst machinery, projecting to first-order relay nuclei whose own A-current retains Kv4.3 backup. Among glutamatergic thalamic cells, the higher-order/midline nuclei (PVT, LP, Xi, intralaminar) are where V404M should bite hardest and the first-order sensory relays least.

## Data and methods

ABC Atlas WMB-10Xv2/v3 (release 20250531; log2 CPM+1 matrices), cells from the TH dissections plus RT-ZI, ZI, LGv/IGL, SPA/PoT and pretectal subclasses and TH/MH-LH classes from the HY and MB dissections (282,078 cells, 183,553 neurons after neurotransmitter filtering). Loop elements defined from CCFv3 nucleus labels (THALMANAC `annotations_c2n_combined.csv`) with subclass fallback (`thalamus/th_common.py`). Because 10x counts include intronic reads, the very long Kv4-complex genes (Kcnip4 1.2 Mb, Dpp10 1.4 Mb, Dpp6 1 Mb, Kcnd2 0.5 Mb) are inflated relative to each other; cross-gene ratios were calibrated with one factor per gene anchoring the whole-thalamus 10x mean to the bulk thalamus WT exon CPM (`thalamus/calibrate_th10x.py`). Within-gene contrasts between cell types are uncalibrated 10x values. Hartley subnetworks: per-cell Spp1 > 0 (415 of 9,075 GABAergic RT-ZI cells; Spp1 mRNA is sparse in 10x) vs rest (`thalamus/trn_subnetworks.py`). Scripts, tables and figures in `thalamus/` and `tables/th10x_*`, `figures/th10x_*` of the Kv4-isoform-survey repository.

## Results

### 1. Kv4 paralogs by loop element (calibrated; Kcnd2 share of Kv4 mRNA)

| loop element | n cells | Kcnd2 share | Kcnd3 raw CPM | Kcnd2⁺ cells | KChIP partners | DPP6 share | KChIP:Kv4 | DPP:Kv4 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TRN (RT) | 9,082 | 0.71 (0.77 RT-annotated) | 39 | 96% | KChIP1 57, KChIP2 36 | 45% | 3.1 | 8.4 |
| first-order sensory relay (VPM/VPL, LGd, MG) | 30,799 | 0.43 | 343 | 100% | KChIP3 61, KChIP4 35 | 72% | 1.0 | 1.6 |
| higher-order relay (PO, LP, MD, SGN) | 24,533 | 0.53 | 214 | 100% | KChIP3 53, KChIP4 38 | 71% | 1.4 | 1.7 |
| motor relay (VAL, VM) | 3,153 | 0.47 | 213 | 100% | KChIP4 75 | 63% | 0.7 | 2.2 |
| intralaminar / midline | 29,838 | 0.57 | 144 | 100% | KChIP4 51, KChIP3 21, KChIP1 18 | 66% | 1.1 | 2.7 |
| anterior / limbic (AD, AV, AM, LD, PVT) | 26,773 | 0.47 | 265 | 100% | KChIP4 54, KChIP1 22 | 69% | 0.9 | 2.1 |
| epithalamus (MH, LH) | 6,223 | 0.47 | 331 | 98% | KChIP1 60 | 69% | 0.7 | 2.1 |
| ZI / LGv / IGL | 37,781 | 0.54 | 119 | 99% | KChIP1 51, KChIP2 19 | 72% | 1.2 | 3.6 |

Per-nucleus values for 50 nuclei with ≥300 cells are in `tables/th10x_calibrated_nuclei.csv`; the Kcnd2-share ranking runs from RT 0.77 and LP/SGN, Xi, PVT/IMD ~0.70 down to LGd 0.31, VPM/VPL 0.30, LD 0.28.

### 2. TRN core vs shell (Hartley et al. 2024)

| | Spp1⁺ core, n = 415 | rest of RT, n = 8,660 |
| --- | --- | --- |
| Kcnd2 / Kcnd3 (raw CPM) | 1,295 / **7** (Kcnd3 detected 20%) | 1,150 / 38 (32%) |
| Kcnip1 / Kcnip2 / Kcnip4 | **830** / 38 / **109** | 410 / 55 / 450 |
| Dpp10 / Dpp6 | 3,020 / 860 | 2,250 / 950 |
| Pvalb / Sst | 2,840 / 613 | 1,610–2,070 / 830–860 |
| Cacna1h / Cacna1i | 47 / 192 (det. 59% / 93%) | 14 / 85 (19% / 63%) |

The differences hold within Pvalb-high cells (Kcnip1 813 vs 498, Kcnip4 118 vs 323, Kcnd3 11 vs 24), so they are properties of the Spp1⁺ cells, not of Pvalb level. Spp1 detection at 4.6% certainly under-calls core cells, so the "rest" group is diluted with undetected core cells and the true contrasts are larger. 340 of the 415 Spp1⁺ cells are in supertype Gnb3_5 (cluster 1604: 79% Spp1⁺).

### 3. Predicted TRN variant composition (proxy from cortical interneurons; measured values pending)

TRN shares lineage (Nkx2.1) and gene-level complex with cortical Pvalb/Sst interneurons, whose variant fractions were measured in Tasic 2016 (full-length): Pvalb cells carry KChIP1c 32% / 1a 23% / 1b 18% / 1e 6% (TM), KChIP2x 78% (TM), DPP10d 86%; Sst cells carry more KChIP1b, KChIP2a and DPP10a. The working prediction for the core TRN complex is therefore Kv4.2 + KChIP1 (myristoyl/palmitoyl mix with some 1e) + KChIP2x + DPP10d. Measurement needs full-length TRN libraries (sorted Viaat-Cre/Spp1-Cre;Ai14 pools, or microdissected TRN bulk); Patch-seq depth is insufficient for these 5′ variants, and 10x/MERFISH cannot see them at all.

### 4. Relay-neuron variants (ThalamoSeq, in progress)

ThalamoSeq (Phillips et al. 2019; 120 pooled full-length libraries from 22 retrogradely labelled cortex-projecting nuclei, single-end 100–105 bp, 30–45 M reads/pool, unstranded) is being quantified with the lab's salmon isoform pipeline. The first pool (PVT) shows reads spread across four Dpp6 and five Kcnip4 N-terminal variants and the Kcnd3 S/L pair, so the 5′ ends are covered and per-nucleus variant fractions with 3–9 replicates per nucleus will follow (`thalamus/thalamoseq_variants.py`). RT is absent from ThalamoSeq (it does not project to cortex).

## Interpretation

1. **TRN is the thalamic cell type without a safety net.** Every relay group keeps 30–70% of its Kv4 pool as Kv4.3, which a heterozygous V404M in Kv4.2 cannot poison; RT keeps 3%. If V404M acts as a dominant-negative within tetramers, RT neurons lose the largest fraction of their A-current of any thalamic cell, and the core cells — the ones that gate first-order sensory and motor relay and carry the most T-type current — lose the most. A-current loss in TRN shortens the interval between bursts and strengthens the rebound, which is the known direction for spike-wave generation.

2. **The two sides of the loop use different complexes.** TRN: Kv4.2–KChIP1/KChIP2–DPP10. Relay: Kv4.2/Kv4.3–KChIP3/KChIP4–DPP6. The dentate/CA1 complex (Kv4.2–KChIP4/KChIP2–DPP6) is a third. The V404M phenotypes therefore involve the mutant subunit in three different auxiliary contexts, and the structural work on how V404M interacts with KChIP and DPP should cover KChIP1/DPP10 (TRN), KChIP3/DPP6 and KChIP4/DPP6 (relay, hippocampus).

3. **Conditional predictions.** Viaat-Cre × cKI puts V404M in RT (93% Slc32a1⁺) and cortical interneurons: spike-wave discharges, TRN and cortical ΔFosB — observed. Vglut2-Cre × cKI puts it in relay neurons and not RT: if relay-cell V404M were sufficient for spike-wave discharges, Vglut2-cKI animals would show them; the fear-conditioning data show no deficit in that line, and its EEG is the direct test. Spp1-Cre × cKI restricts V404M to the core TRN subnetwork and is the sharpest test of the TRN hypothesis; Ecel1-Cre × cKI the complementary one.

## Reliability

| claim | rating | basis |
| --- | --- | --- |
| RT has 10–30× less Kcnd3 than relay neurons | solid | within-gene 10x contrast, 9,082 cells, no calibration needed |
| RT uses KChIP1/KChIP2 + DPP10; relay uses KChIP3/4 + DPP6 | solid | within-gene contrasts; consistent with Serôdio & Rudy 1998 ISH for Kv4.3 vs Kv4.2 |
| Calibrated cross-gene shares (Kcnd2 0.77 in RT etc.) | moderately solid | one calibration factor per gene from bulk exon data; assumes neuron-dominated bulk |
| Spp1⁺ core = purest Kv4.2-only, KChIP1-rich, T-type-rich TRN cells | moderately solid | 415 cells; holds within Pvalb-high; Spp1 under-detection dilutes the contrast |
| FO vs HO relay Kcnd2 gradient | solid | per-nucleus, thousands of cells each |
| TRN variant composition | prediction only | cortical interneuron proxy; needs full-length TRN libraries |

## Proposed experiments

1. **EEG in Vglut2-Cre × cKI** (animals exist): relay-cell V404M alone — spike-wave or not. Together with the Viaat-Cre result this assigns the loop phenotype to one side.
2. **Spp1-Cre and Ecel1-Cre × cKI** (Hartley lines): core vs shell TRN; EEG and TRN ΔFosB by sector.
3. **Full-length TRN RNA-seq**: sorted tdTomato⁺ pools from Viaat-Cre;Ai14 and Spp1-Cre;Ai14 (SMART-seq v4, 20–40 M reads, 3–4 pools per line), giving the measured KChIP1/KChIP2/DPP10 variant fractions and, from the same cells, the cKI recombination efficiency.
4. **TRN A-current recordings in V404M/+ vs +/+ slices**, core (Pvalb-high, Spp1-reporter) cells first: the direct test of the safety-net argument.
5. **Structural/MD**: V404M in Kv4.2 complexes with KChIP1 + DPP10 alongside the KChIP3/4 + DPP6 complexes already modelled.

## Provenance

`thalamus/extract_th_10x.py` → `summarize_th_10x.py` → `calibrate_th10x.py`, `trn_subnetworks.py` (run on PFAFF-LAMBDA, cache E:\Claude_workspace\abc_atlas); tables `th10x_loop_summary.csv`, `th10x_nucleus_summary.csv`, `th10x_cluster_means.csv`, `th10x_calibrated_*.csv`, `th10x_trn_core_vs_shell.csv`; figures `th10x_kv4_calibrated.png`, `th10x_kv4_loop_heatmap.png`. ThalamoSeq: `pipeline/build_thalamoseq_runlist.py`, `run_thalamoseq_isoforms.sh`, `collapse_thalamoseq.py`; `thalamus/thalamoseq_variants.py`. THALMANAC: AllenNeuralDynamics/thalmanac (abc-merfish-analysis). Literature from memory, verify before quoting: Hartley et al. 2024 Cell Rep 43:115037 (Spp1/Ecel1 TRN subnetworks); Phillips et al. 2019 Nat Neurosci 22:1925 (ThalamoSeq); Serôdio & Rudy 1998 J Neurophysiol (Kv4 subunit ISH); Clemente-Perez et al. 2017 Cell Rep (PV vs SOM TRN cells and T-type currents).
