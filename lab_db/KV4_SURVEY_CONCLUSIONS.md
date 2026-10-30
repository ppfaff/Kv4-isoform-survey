# Kv4 Auxiliary-Subunit Isoform Survey: Conclusions and Reliability

Oct 29, 2026 · @Paul

> **Revision note (Oct 30, 2026).** The dendrite-vs-soma compartment analysis (companion document *Kv4 complex mRNAs are somatic*, Farris 2019 laminar LCM + Middleton 2019 single-neuron data) supersedes one interpretation in this document. The somatic Kv4-null cells (DG 18%, Vip 27%, Pvalb 20%) were described below as "post-transcriptional, most plausibly dendritic localization of the mRNA". In adult hippocampus Kcnd2, Kcnip4, Dpp6 and Kcnip2 mRNAs are depleted 2–4× in the dendritic layers while canonical dendritic mRNAs are enriched 3–9×, and the Kcnd2 signal in the neuropil is 52% intronic (nuclear pre-mRNA from resident interneurons and glia). These cells therefore have **low Kv4 mRNA, not relocated Kv4 mRNA**, and should be expected to have smaller A-currents; the Kv4-null group shows no Kcnd1/Kcnd3 compensation. The alternatives still open are splicing efficiency, nuclear export and mRNA stability. Two further points from that analysis bear on this document: salmon underestimates Kv4.3L ~4-fold in bulk data (junction-level exon inclusion 0.17 vs salmon 0.04 in cell-body layers), so the Kv4.3L fractions here are lower bounds, consistent with the "flagged" rating already given to the Kv4.3S/3L balance; and the Kcnip3 cys-rich variant behaves like a nuclear-retained species (28% in SSv4 cells → 19% in cell-body layers → 4% in neuropil), so its functional status is doubtful. Everything else in this document stands.

## Summary

Across ~11,300 single neurons from mouse V1 and hippocampal formation, the Kv4 channel complex (Kcnd1–3 pore subunits, Kcnip1–4 KChIPs, Dpp6/Dpp10) is assembled from cell-class-specific splice variants, and the pore subunit — not the auxiliary subunits — is the variable, regulated component. Interneurons and principal cells use nearly disjoint KChIP programs (KChIP1 vs KChIP2–4), a split that holds in both regions; the DPP split is by variant and region instead — DPP10d plus DPP6 in interneurons, DPP10a plus DPP6 in cortical pyramidal cells, and DPP6 alone in hippocampal principal cells. KChIPs and DPPs are expressed in almost every neuron; Kv4 mRNA is absent from the soma of a type-specific fraction of cells (DG granule cells 18%, Vip 27%, Pvalb 20%), and the Allen intron/exon data show those cells are still transcribing Kcnd2 — the absence is post-transcriptional, most plausibly dendritic localization of the mRNA. *[Revised Oct 30: not dendritic localization — see revision note at top.]*

Three findings are new to this survey: hippocampus-specific Sst and Lamp5 Lhx6 interneuron types that have switched the KChIP genes off; near-projecting (NP) pyramidal cells in both cortex and subiculum that run an interneuron-style KChIP1 complement on a ~92% Kv4.2 pore; and a set of regional variant differences at matched chemistry (DPP10a is neocortical; hippocampal interneurons shift toward Kv4.3S and KChIP1b). A fourth result is methodological and matters for everything else: the two V1 datasets differ in read chemistry, and that alone shifts the apparent fractions of several variants (Kv4.3L, DPP6L, KChIP1a/b, KChIP2a/x, DPP10c) by 4–14 points in every cell class, so cross-dataset isoform fractions are not comparable without a same-chemistry control. The control run (1,511 Allen VISp cells) was completed on 2026-10-29 and underlies every regional claim below.

For Kv4.2-V404M, the cells most exposed are near-projecting L5/6 and subicular pyramidal cells, CA1/CA3/subicular pyramidal cells, Sst and Pvalb interneurons, and DG granule cells; the least exposed are Vip, Sncg and Lamp5 interneurons, which build their A-current on Kv4.3.

## Datasets and what was compared

Three SMART-seq datasets were quantified with one pipeline (salmon 1.x selective alignment against a GENCODE vM36 transcriptome-only index; per-transcript reads and TPM for the 9 genes; variants assigned by a curated transcript-to-isoform map), then compared at the level of variant fractions, gene and family expression, and per-cell completeness of the Kv4–KChIP–DPP complex.

| Dataset | Region | Cells quantified | Chemistry | Depth (median reads/cell) | Role in this survey |
| --- | --- | --- | --- | --- | --- |
| Tasic 2016 (GSE71585) | V1 (VISp) | 1,677 of 1,679 annotated; 49 types | SMART-Seq v1, 100 bp single-end | 8.7 M | Pilot dataset; all analyses developed here (first on a 460-cell preview, then the full set) |
| Allen SSv4 hippocampal formation (Yao 2021) | HIP, CA fields, DG, SUB-ProS, adjoining cortex | 8,131 of 8,137 (26 subclasses; 18 with ≥15 cells) | SMART-Seq v4, 50 bp paired-end | 0.70 M | Main regional comparison; 12× shallower than Tasic 2016 |
| Allen SSv4 VISp (Tasic 2018 cells) | V1 | 1,511 of a stratified subset (100 per subclass) | SMART-Seq v4, 50 bp paired-end | ~2.5 M | Chemistry control: same region as Tasic 2016, same chemistry as the hippocampal set |

The comparisons made, in order:

1. Variant composition per cell type (pooled-read fractions of each family) and per-gene variant fractions, with Kruskal–Wallis / Mann–Whitney tests across types and classes.
2. Per-cell expression of each family (reads per million processed) and of each gene, cells clustered by type, with explicit marking of cells with zero reads.
3. An “A-current complex index” per cell (each family’s CPM divided by its mean over all neurons, then multiplied), and the fraction of depleted cells per type, attributed to the limiting family.
4. Per-cluster breakdown of depleted cells within the variable subclasses (Lamp5, Sst, Vip, Pvalb, Sncg, DG, CT SUB).
5. Cross-check of the salmon gene-level counts against the Allen Institute’s released RSEM exon counts for the same cells (Spearman ρ 0.97–0.999 for 8 of 9 genes; Dpp10 0.87), and use of the Allen combined (exon+intron) matrix to separate nascent transcription from spliced mRNA.
6. Maturity markers (Dcx, Calb2, Sox11, Calb1, Prox1) in DG granule cells, from the Allen matrix.
7. Region comparison in matched classes (Lamp5, Sncg, Vip, Sst, Pvalb, L2/3 IT, L6 CT, L6b), with V1 reads binomially thinned to hippocampal depth for zero-rate comparisons, and the three-way Tasic-2016 / Allen-VISp / Allen-HPF comparison to separate chemistry from biology.
8. Driver-gene expression (Tshz2, Slc17a8, Slc17a6, Camk2a, Rbp4, Trhr, Fezf2 and others) in NP cells versus all other classes, from the Allen exon matrix.

## How reliability is rated

Each conclusion below carries one of three ratings, decided by whether it survives the specific failure modes found in this project.

- **Reliable** — seen within a single dataset (so chemistry cannot explain it), reproduced in at least one other dataset or region, and robust to depth (holds at full depth and after thinning, or does not depend on zero counts).
- **Tentative** — supported by the data but resting on one dataset, on small cell numbers (<30), on a cross-dataset comparison of isoform fractions, or on an interpretation (for example, mRNA localization) that the data constrain but do not prove.
- **Flagged** — known or likely to be an artifact, or contradicted by the control; kept in the record so the pilot analyses can be re-read correctly.

Two project-specific rules follow from the cross-checks. First, any isoform fraction compared *between* Tasic 2016 and either Allen dataset is at best tentative, because the read chemistry alone shifts several variants by 4–14 points. Second, any statement about cells “lacking” a gene refers to spliced mRNA in the soma of a dissociated cell; for Kcnd2 in particular that is not the same as the gene being off.

## Conclusions, with reliability

### Variant usage by cell class

| Conclusion | Evidence | Rating |
| --- | --- | --- |
| Interneurons and principal cells use nearly disjoint KChIP programs: interneurons KChIP1 (1b > 1a), principal cells KChIP2/3/4, with KChIP1 absent from ~90% of principal cells (NP cells excepted) | Same split in Tasic 2016 V1, Allen VISp and Allen HPF; gene level, so chemistry-independent | Reliable |
| DPP usage splits by variant and region, not by class: interneurons co-express DPP10d and DPP6 (DPP10 dominant only in hippocampal Lamp5, 64%); cortical pyramidal cells are DPP10a-dominant (L2/3 IT 0.39–0.46, L5 IT and L6 CT 0.7–0.9); hippocampal CA/DG/SUB cells are nearly pure DPP6 (DPP10 <5%) | All three datasets; the DPP10a-vs-DPP6 contrast between cortical and hippocampal principal cells is gene-level and chemistry-independent | Reliable |
| Lamp5, Sncg and Vip interneurons are Kv4.3S-dominated (49–69% of Kv4); Sst and Pvalb are Kv4.2-dominated (50–82%) | All three datasets | Reliable |
| Hippocampal principal cells are strongly Kv4.2: CA1 84%, SUB 79%, NP SUB 92%, DG 71%, CA3 56% | Allen HPF, 8,131 cells | Reliable |
| KChIP4bL is >90% of all Kcnip4 in hippocampal principal cells, versus 20–40% in Vip/Sst | Allen HPF | Reliable (within dataset) |
| DPP6a is a fixed ~31–40% of Dpp6 in every hippocampal subclass, suggesting an unregulated splice ratio | Allen HPF, all 18 subclasses | Tentative (one dataset; the Allen chemistry assigns ~9 points more of Dpp6 reads to DPP6L in every class, which dilutes the DPP6a fraction by a constant amount and could itself make the ratio look fixed) |
| Kv4.3L share depends on both normalization and chemistry: by pooled reads it is 2–10% of Kv4 in HPF, but by TPM (molecule-weighted) 13–30% in Lamp5/Vip/Sst/DG and 4–8% in CA1/NP SUB/Pvalb, because the Kv4.3L transcript model lacks the 7-kb 3′UTR and loses the shared-UTR reads to Kv4.3S; in Tasic 2016 it is 0–3% by either measure | Allen VISp and HPF agree with each other; Tasic 2016 disagrees with both; the TPM values are the molecule-level estimate | Flagged: Kv4.2 and Kv4.1 shares are unaffected by the weighting (≤2 points), so the Kv4.3S/3L balance is the only Kv4 number that is not settled |
| Vip is the one interneuron class built on KChIP4 (50–66% of its KChIP); Pvalb is the one with substantial KChIP2 | Allen HPF and Allen VISp agree | Reliable |

![Figure 1. Hippocampal variant composition by subclass](https://github.com/ppfaff/Kv4-isoform-survey/blob/main/figures/hip_variant_scheme_by_subclass.png?raw=true)

**Figure 1. Variant composition per hippocampal subclass (Allen SSv4, 8,131 cells).** Each row is a subclass; the three panels are the Kv4, KChIP and DPP families, each bar the pooled-read fraction of every variant in that family, drawn most-common segment first. Rows are grouped inhibitory (purple labels) / excitatory (blue) / non-neuronal (grey); the right column gives cells per subclass. Methods: salmon quantification against a GENCODE vM36 transcriptome-only index; variants assigned by the transcript map in the Appendix; retained-intron transcripts dropped only when unassignable; subclasses with fewer than 15 cells omitted; a grey bar marks fewer than 20 reads on the family. Fractions are read-weighted and therefore length-biased for Kv4.3S vs 3L (see Appendix); all other variant pairs are within 2 points of their TPM-weighted values.

![Figure 2. V1 variant composition by type, Tasic 2016](https://github.com/ppfaff/Kv4-isoform-survey/blob/main/figures/ctx_variant_scheme_by_type.png?raw=true)

**Figure 2. The same analysis for the 49 Tasic 2016 V1 types (1,677 cells, 100 bp single-end).** Layout and methods as Figure 1 with a 5-cell minimum per type. Because of the chemistry effects documented in the section on the two V1 datasets, minor-variant percentages here are not directly comparable with Figure 1; the dominant variant per type is.

### Family-level expression and complete complexes

| Conclusion | Evidence | Rating |
| --- | --- | --- |
| KChIPs and DPPs are expressed in nearly every neuron (family null rates ≤2% in most subclasses) and their levels do not track Kv4 (Spearman ρ ≈ 0.07–0.09) | Allen HPF; same pattern in the V1 pilot | Reliable |
| CA1, CA2, CA3, SUB-ProS and NP SUB pyramidal cells carry a complete Kv4–KChIP–DPP complement with the least cell-to-cell variability (6–13% depleted by the index) | Allen HPF | Reliable |
| Kv4 is the variable, type-specific component: somatic Kv4-null cells are 18% of DG, 27% of Vip, 20% of Pvalb, 12% of Lamp5, 10% of Sst, ≤2% of CA/SUB | Allen HPF; V1 matched classes give the same ranking after depth thinning (Vip 38%, Pvalb 26%, Sst 10%, L2/3 4%) | Reliable |
| Those null fractions are not a depth or quality effect | No depth difference between null and expressing cells (Mann–Whitney p = 0.98 GABA, 0.08 Glut); V1 Vip cells are still 25% null at 8.7 M reads | Reliable |

![Figure 3. Per-cell Kv4 signal clustered by subclass](https://github.com/ppfaff/Kv4-isoform-survey/blob/main/figures/hip_percell_kv4_by_subclass.png?raw=true)

**Figure 3. Kv4 pore-subunit mRNA per cell, hippocampal neurons clustered by subclass.** Each column is one cell; coloured points are that cell’s Kv4.1, Kv4.2, Kv4.3S and Kv4.3L reads per million reads processed (CPM, log scale; zeros not plotted). Subclass blocks are ordered left to right by the subclass mean (black step), and cells within a block by their own total. Gray crosses at the bottom mark cells with zero reads on all three Kcnd genes (“somatic Kv4-null”). The colour strip gives subclass identity. Methods: CPM uses reads processed by salmon as the per-cell denominator; non-neuronal cells (37, of which 36 are Kv4-null) are omitted. The null fraction is the quantity carried through the rest of the document; its interpretation is addressed in Figure 8.

![Figure 4. Per-cell KChIP gene expression](https://github.com/ppfaff/Kv4-isoform-survey/blob/main/figures/hip_percell_KChIP_genes_by_subclass.png?raw=true)

**Figure 4. KChIP family at gene level per cell (Kcnip1–4; all transcripts of each gene).** Same layout as Figure 3 with blocks ordered by subclass-mean KChIP. The near-absence of KChIP-null cells in principal cells, the Sst and Lamp5 blocks sitting lower than every other interneuron block, and the class-wide KChIP1 (interneuron) versus KChIP2/3/4 (principal cell) switch are read directly from this figure.

![Figure 5. Per-cell DPP gene expression](https://github.com/ppfaff/Kv4-isoform-survey/blob/main/figures/hip_percell_DPP_genes_by_subclass.png?raw=true)

**Figure 5. DPP family at gene level per cell (Dpp6, Dpp10).** As Figure 4. Principal cells form a tight DPP6 band with DPP10 scattered far below and absent in roughly half the cells; interneurons co-express both.

![Figure 6. A-current complex index by subclass](https://github.com/ppfaff/Kv4-isoform-survey/blob/main/figures/hip_acurrent_index_by_subclass.png?raw=true)

**Figure 6. A-current complex index.** Left: per-subclass distribution (10–90% whiskers, 25–75% bar, median) of the index = (Kv4 CPM / mean) × (KChIP CPM / mean) × (DPP CPM / mean), with each family’s mean taken over all 8,093 hippocampal neurons (Kv4 86, KChIP 333, DPP 127 CPM), so a cell at the neuron average for all three scores 1; the grey number is cells at zero. Right: percentage of cells with index < 0.1, stacked by which family is lowest in that cell. Methods note: the 0.1 cutoff is a convention and the three families have different per-cell variance, so the percentages are for ranking subclasses, not absolute prevalence; the fractions at exactly zero and at < 0.01 are in `hip_acurrent_index_by_subclass.csv`.

![Figure 7. Kv4, KChIP and DPP totals per cell ordered by Kv4](https://github.com/ppfaff/Kv4-isoform-survey/blob/main/figures/hip_kv4_kchip_dpp_total.png?raw=true)

**Figure 7. Family totals per cell with cells ordered by Kv4.** Black = Kv4, green = KChIP, red = DPP (CPM, log scale). KChIP and DPP stay flat across the entire Kv4 gradient and past the point where Kv4 reaches zero (Spearman ρ with Kv4 = 0.07–0.09 in both classes). This is the basis for the statement that the auxiliary subunits are near-constitutive while the pore subunit varies.

### The intron/exon finding on Kv4-null cells

| Conclusion | Evidence | Rating |
| --- | --- | --- |
| Our salmon quantification agrees with the Allen RSEM exon counts cell by cell (Kcnd2 ρ = 0.997, medians 22 vs 21 reads; 434 of 446 DG null cells are null in both) | Allen exon-only matrix | Reliable (pipeline validated) |
| Kv4-null cells are still transcribing Kcnd2: intronic (pre-mRNA) Kcnd2 reads are ~400 CPM in null, low and high DG cells alike, and exonic vs intronic Kcnd2 are uncorrelated (ρ = 0.01) | Allen combined vs exon matrices | Reliable |
| In every subclass, most somatic Kv4-null cells still carry intronic Kv4 reads (nuclear pre-mRNA, i.e. the gene is being transcribed); cells with neither spliced nor intronic reads are ≤6% (Vip) and ~0% in principal cells. Gene-level counts cannot say whether a particular intron is retained; that needs read-level coverage of the Kcnd2 locus | Allen HPF, all subclasses | Reliable |
| The effect is specific to Kcnd2: the same cells have normal exonic and intronic Kcnip4 and Dpp6, normal exonic yield, and only a marginally higher nuclear read fraction | Allen HPF DG | Reliable |
| Explanation: spliced Kcnd2 mRNA is dendritically localized and lost with the dendrites in dissociation, so “null” means “none in the soma”; the continuous exonic spread at constant intronic level reads as a per-cell somatic-retention phenotype | Consistent with the known dendritic targeting of Kv4.2 mRNA, but not demonstrated here. Competing explanations: rapid somatic decay (indistinguishable in this data); regulated nuclear retention (argued against: null cells do not have a larger intronic pool); failed elongation of the 480-kb gene; or intronic reads from a nested transcript rather than pre-mRNA (which would also explain the zero exon–intron correlation and the complete absence of exonic reads in null cells). Per-intron coverage (item 5b′) and exon/intron smFISH (item 1) separate these | Tentative |
| The DG null fraction is not developmental: no Dcx+/Calb2+/Calb1− cells exist in the set; null cells are Calb1+ (98%), Prox1+ (96%); no marker correlates with Kcnd2 | Allen exon matrix | Reliable (for this adult dataset) |

![Figure 8. DG granule cells: exonic vs intronic Kcnd2](https://github.com/ppfaff/Kv4-isoform-survey/blob/main/figures/dg_kcnd2_exon_vs_intron.png?raw=true)

**Figure 8. Kv4-null DG granule cells transcribe Kcnd2 and are mature.** Left three panels: per DG cell (n = 2,458), exonic (spliced mRNA) versus intronic (pre-mRNA) reads for Kcnd2, Kcnip4 and Dpp6, CPM on log axes; cells at exonic zero are drawn at the left edge. Colour is the cell’s Kv4 group from our salmon quantification: null = zero reads on Kcnd1/2/3; low/high = below/above the median Kcnd2 among expressing cells. Right: percentage of cells in each group with each maturity or identity marker detected. Methods: exonic counts are the Allen Institute’s released RSEM exon-only matrix for the same cells; intronic = the Allen combined (exon+intron) matrix minus the exon-only matrix, clipped at zero; both divided by reads processed. Markers come from the same matrices. The 446 null cells have intronic Kcnd2 of ~400 CPM, identical to the other groups (p = 0.6), while Kcnip4 and Dpp6 are unchanged in both compartments; exonic and intronic Kcnd2 are uncorrelated (ρ = 0.01). No Dcx+/Calb2+/Calb1− cells exist in this adult dataset. Caveat: gene-level intron counts cannot localize the signal within the 480-kb gene; see the alternatives in the table above.

### KChIP-off interneuron types

| Conclusion | Evidence | Rating |
| --- | --- | --- |
| Hippocampus-specific Sst types (77, 78, 102–105 Sst HPF) have little or no KChIP; 77\_Sst HPF has a median of zero with 82% of cells lacking nascent KChIP transcription | Allen HPF; nascent and spliced agree | Reliable |
| Two Lhx6+ Lamp5 clusters (7\_ and 8\_Lamp5 Lhx6) are KChIP-low to KChIP-off (medians 68 and 36 CPM vs 629 and 416 in 5\_/6\_Lamp5 Lhx6), transcriptionally | Allen HPF; 83 of 92 KChIP-null Lamp5 cells are in these two clusters; intronic Kcnip1 median 0 | Reliable |
| In Lamp5 the Kv4 drop-out and the KChIP drop-out are independent (9 cells with both vs 10.5 expected) | Allen HPF | Reliable |
| Sncg and Lamp5 share a variant recipe but differ in reliability: Sncg types are 0–17% depleted, Lamp5 types 23–70%, with the Lhx6+ (MGE-derived) Lamp5 types the worst | Allen HPF per cluster | Reliable |
| Cortical Smad3 cells correspond to the Lhx6+ Lamp5 type by their KChIP4-shifted, KChIP1-poor profile | Tasic 2016 (13 cells) matched by profile | Tentative |

![Figure 9. Lamp5 interneurons: family totals by Kv4 order and by cluster](https://github.com/ppfaff/Kv4-isoform-survey/blob/main/figures/hip_lamp5_family_totals.png?raw=true)

**Figure 9. Lamp5 interneurons (n = 891): Kv4, KChIP and DPP per cell.** A: cells ordered by total Kv4; the red dashed line marks the start of the zero-Kv4 cells; crosses along the bottom mark cells with zero reads for each family. KChIP and DPP are flat across the Kv4 decline and past it. B: the same cells grouped by Allen cluster (clusters ordered by median KChIP; cells within a cluster by Kv4), with thick bars at cluster medians (green KChIP, black Kv4, red DPP). The KChIP deficit is a cluster property: 7\_ and 8\_Lamp5 Lhx6 have medians of 68 and 36 CPM against 629 and 416 in 5\_ and 6\_Lamp5 Lhx6, and carry 83 of the 92 KChIP-null Lamp5 cells; the intron data (not shown) give zero nascent Kcnip1 in those clusters, so the silencing is transcriptional. Methods: CPM by reads processed; cluster labels are Allen’s. The same per-cluster analysis for Sst (table in `hip_exon_vs_intron_by_subclass.csv` and the per-cluster output) identifies the Sst HPF clusters as the KChIP-off Sst types.

### Near-projecting (NP) pyramidal cells

| Conclusion | Evidence | Rating |
| --- | --- | --- |
| NP cells in cortex and subiculum are the only pyramidal class with an interneuron-style KChIP1 complement: KChIP1 is >25% of KChIP in 81% (VISp NP) and 87% (NP SUB) of cells, versus 0–2.5% in every other pyramidal class | Allen VISp and Allen HPF, same chemistry, independently sampled | Reliable |
| NP cells are the most Kv4.2-dominant class (~90–92% of Kv4) with the highest Kv4 expression, highest complex index (median 2.3) and lowest depleted fraction (6%) among hippocampal types | Allen HPF and Allen VISp | Reliable |
| The Tasic 2016 “L5 Ucma” type is not the NP population (zero KChIP1; a KChIP3a/DPP10a profile like L5 IT) | Expression profile | Reliable (as a negative) |
| Driver coverage: Tshz2 in 100% of NP cells (also SUB-ProS, CT SUB); Slc17a8 in 67% of cortical NP but 0% of NP SUB; Camk2a detected in 71–78% but at the lowest level of any pyramidal class; Slc17a6 in 83% of NP SUB, 44% (low) of cortical NP; Rbp4 in 78% of cortical NP | Allen exon matrix | Reliable (expression); line behavior is inferred |

### Regional differences at matched chemistry (Allen VISp vs Allen HPF)

| Conclusion | Evidence | Rating |
| --- | --- | --- |
| DPP10a is a neocortical pyramidal-cell variant and nearly absent from hippocampal-formation L2/3 IT (0.39 vs 0.06), replaced by DPP6S and DPP6L | BH-corrected p < 0.001 | Reliable (note the HPF L2/3 cells are retrohippocampal cortex, so this is cortex-vs-cortex) |
| Hippocampal Lamp5 and Sncg use less KChIP4bL (0.30 → 0.09; 0.26 → 0.05) and more KChIP1b; Vip switches KChIP4 variant (KChIP4d 0.09 → 0.30) | BH p < 0.001 | Reliable |
| Hippocampal interneurons shift toward Kv4.3S: Lamp5 0.38 → 0.51, Vip 0.48 → 0.67, Pvalb 0.10 → 0.27 | Same direction in all three classes | Reliable |
| Sst composition is nearly identical between regions; the Sst regional story is KChIP amount, not variants | Allen VISp vs HPF | Reliable |
| Hippocampal interneurons have ~2× the KChIP expression of cortical ones (Vip 452 vs 203 CPM, Pvalb 359 vs 200) | Allen HPF vs Tasic 2016 CPM | Tentative (cross-chemistry; to be confirmed with Allen VISp depth data) |
| L6 CT and L6b regional differences | 17–30 hippocampal cells | Tentative (small n) |

![Figure 10. Three-way variant comparison: Tasic 2016 V1, Allen VISp, Allen HPF](https://github.com/ppfaff/Kv4-isoform-survey/blob/main/figures/threeway_variant_composition.png?raw=true)

**Figure 10. Variant composition in matched cell classes across three datasets.** Each class has up to three bars: grey label = Tasic 2016 VISp (100 bp single-end), blue = Allen SSv4 VISp (50 bp paired-end, 100 cells per subclass), red = Allen SSv4 hippocampal formation. Grey-vs-blue isolates the chemistry effect (same region, same types); blue-vs-red is the regional comparison at matched chemistry. Classes below the divider (L4, L5 IT, L5 PT, L5/6 NP, L6 IT) have no hippocampal counterpart and show the chemistry comparison only. Methods: Tasic 2016 types were mapped to Allen subclasses as listed in the two-V1-datasets section; pooled-read fractions as in Figure 1; the significance tests in `threeway_tests.csv` use per-cell fractions in cells with ≥20 reads on the family (Mann–Whitney, Benjamini–Hochberg across all tests). The L5/6 NP row’s grey bar (Tasic “L5 Ucma”) is not the NP population and should be disregarded; the blue NP bar, with its KChIP1-based complement, stands on its own and matches hippocampal NP SUB. The regional-difference table above and the chemistry-shift table in the two-V1-datasets section are both read from this figure.

## The two V1 datasets: consistent and inconsistent

Tasic 2016 (100 bp single-end, ~8.7 M reads/cell) and the Allen SSv4 VISp cells (50 bp paired-end, ~2.5 M reads/cell) sample the same region and the same cell types. Comparing 13 matched classes (8 shared with hippocampus plus L4, L5 IT, L5 PT, L5/6 NP, L6 IT) isolates what the chemistry does to the isoform assignment. Figure 10 shows all three datasets side by side for each class.

### Consistent between the two V1 datasets

- **Which variant dominates in each class.** The Kv4.2-vs-Kv4.3S assignment of every interneuron and pyramidal class, the KChIP1 program of interneurons versus KChIP2–4 in principal cells, the DPP10d-plus-DPP6 program of interneurons versus DPP10a-plus-DPP6 in pyramidal cells, and type-level signatures such as KChIP3a in L5 IT and L6 CT, KChIP4bL+KChIP3a in L5 PT, DPP10a in L2/3–L5 pyramidal cells. Examples of agreement: Lamp5 Kv4.2 0.53 vs 0.50; Vip Kv4.3S 0.47 vs 0.48; L2/3 IT DPP10a 0.46 vs 0.39.
- **The somatic Kv4-null phenomenon and its ranking across classes.** Vip and Pvalb high, Sst and L2/3 low, in both; and it persists at 8.7 M reads (25% of Tasic Vip cells), so it is not a depth artifact.
- **Gene-level expression patterns** (which KChIP and DPP genes each class expresses), including the absence of KChIP1 from pyramidal classes other than NP.

### Inconsistent between the two V1 datasets (chemistry effects)

The Allen 50-bp paired-end chemistry, relative to Tasic 2016, shifts the following fractions with the same sign in essentially every class (mean shift over 13 classes):

| Variant | Shift (Allen − Tasic 2016) | Classes with the same sign | Likely cause |
| --- | --- | --- | --- |
| KChIP1b | +11 points | 10 of 13 | KChIP1a/1b differ by a 33-nt insert; read length and pairing change how spanning reads are assigned |
| DPP6L | +9 | 13 of 13 | Alternative 5′ exon; DPP6L is near zero in Tasic 2016 |
| KChIP1a | +7 | 13 of 13 | As for KChIP1b |
| KChIP2a | +5 | 13 of 13 | KChIP2 N-terminal exons |
| Kv4.3L | +4 | 13 of 13 | 57-nt C-terminal insert; Kv4.3L is 0–3% in Tasic 2016, 7–19% in Allen data |
| KChIP3a | −11 | 13 of 13 | Reciprocal of the KChIP1/2 gains in classes that also express them |
| Kv4.2 | −7 | 11 of 13 | Reciprocal of the Kv4.3L/3S gains |
| DPP10c | −5 | 13 of 13 | Dpp10-202 carries a 27-kb annotated transcript; its UTR collects reads differently by chemistry |
| KChIP2x, KChIP4d, KChIP1c | −3 to −5 | 8–12 of 13 | Minor variants distinguished by short segments |

Two consequences. Every “regional difference” in Kv4.3L, DPP6L, DPP10c or KChIP2x first seen between Tasic 2016 and the hippocampal data was chemistry, not biology, and is withdrawn. And the pilot figures built on Tasic 2016 remain correct about which variants each type uses, but their minor-variant percentages should not be quoted as absolute values; the Allen-chemistry values are the better estimate where the two disagree, because paired-end reads resolve short alternative segments more reliably.

One mapping inconsistency is not chemistry: the Tasic 2016 “L5 Ucma” type, mapped to L5/6 NP, has none of the NP signature (zero KChIP1, KChIP3a/DPP10a profile). It is an IT-like population, and the NP class is effectively absent from the 2016 annotation used here.

## Cells most and least sensitive to Kv4.2-V404M

Exposure to a Kcnd2 mutation is taken here as the product of three things the data measure: how much of a cell’s A-current pore is Kv4.2 rather than Kv4.1/4.3, how much Kv4 it expresses overall, and how consistently the cell carries a complete complex. Where it is known, the auxiliary-subunit context is noted, because KChIP and DPP6 binding strongly modulate the inactivation gating that V404M disrupts. The Kv4-null caveat applies throughout: cells with no somatic Kcnd2 mRNA are still transcribing it, so they are not spared. The Kv4.2 shares come from Figures 1 and 10, the expression levels and null fractions from Figure 3, and the completeness measures from Figure 6.

### Most sensitive

| Cell class | Kv4.2 share of Kv4 | Kv4 level and completeness | Partners | Confidence |
| --- | --- | --- | --- | --- |
| Near-projecting pyramidal cells, L5/6 NP (cortex) and NP SUB (subiculum) | ~90–92% | Highest Kv4 and highest complex index of any hippocampal subclass; 6% depleted; 2% null | KChIP1 (~50%) + KChIP4/3; DPP6 — a unique, interneuron-like complex on a Kv4.2 pore | High: same result in two regions, two independent samples |
| CA1 pyramidal cells (CA1-ProS) | 84% | 113 CPM median; 13% depleted; 1% null | KChIP2a + KChIP4bL + KChIP3a; DPP6S/6a | High |
| Subiculum (SUB-ProS) and CA3 | 79% / 56% | 10% and 9% depleted; ≤2% null | SUB: KChIP4bL 64%; CA3: KChIP2a 44% | High |
| DG granule cells | 71% | 40% depleted, 18% somatic-null, but all transcribing Kcnd2; Kv4.2 in dendrites is the expected site of action | KChIP2a/4bL/3a; DPP6S/6a | High for exposure; the null-cell interpretation is tentative |
| Pvalb interneurons | 64% (HPF), 76–82% (V1) | 32% depleted, 20% null; lowest Kv4 among interneurons after Vip | KChIP1 + KChIP2 (the only interneuron class with substantial KChIP2); DPP10d/DPP6 | Medium-high: Kv4.2-dominant but with a large null fraction |
| Sst interneurons | 50% (HPF), 70% (V1) | 10% null; 49% depleted in HPF, mostly KChIP-limited | KChIP1; in HPF-specific Sst types KChIP is absent, which would leave Kv4.2 channels without KChIP modulation — V404M effects on inactivation could differ there | Medium |
| L5 PT, L5 IT, L4, L6 CT cortical pyramidal cells | 55–75% | Moderate; 17–34% depleted | KChIP3a/4bL; DPP10a | Medium (Kv4.2 share is chemistry-sensitive by ~7 points) |

### Least sensitive

| Cell class | Kv4.2 share of Kv4 | Why least exposed | Confidence |
| --- | --- | --- | --- |
| Vip interneurons | 18% (HPF), 28–44% (V1) | Kv4.3S-based A-current; 27% somatic-null; the lowest Kv4 of any neuronal class; KChIP4-based complex | High |
| Sncg interneurons | 25% | Kv4.3S 65%; very consistent complex (17% depleted) but Kv4.2 is a minority | High |
| Lamp5 interneurons, especially the Lhx6+ types | 32–36% (HPF), ~50% (V1) | Kv4.3S-dominant in hippocampus; the KChIP-off Lhx6+ types carry little functional complex at all | High for hippocampus; in cortex Lamp5 is closer to 50/50 |
| CT SUB | 55% | Small, but 55% depleted and the KChIP3a/DPP10a complex is unusual; exposure is real but the cell population is small | Low (n = 123) |
| Non-neuronal cells | — | No Kv4 expression | High |

Two points for experimental design follow. The classes at the two extremes within one tissue — NP and CA1 pyramidal cells versus Vip and Sncg interneurons — give an internal control: a V404M phenotype that appears in the former and not the latter is attributable to Kv4.2, whereas one that appears in both is not. And because every DG granule cell transcribes Kcnd2 while a fifth keep no spliced message in the soma, DG is the population where dendritic versus somatic effects of the mutation could be distinguished, which smFISH with exon and intron probes would set up.

## Status of the analyses

The analyses were developed on the 460-cell Tasic 2016 preview, then run on the full 1,677-cell set and on the hippocampal data; the Allen VISp control was added last. Not everything has been carried through to every dataset.

| Analysis | Tasic 2016 V1 (full) | Allen HPF | Allen VISp control |
| --- | --- | --- | --- |
| Variant composition by type (v17 scheme figure, tables) | Done (49 types) | Done (18 subclasses) | Done (matched classes only; the full 25-subclass figure not yet made) |
| Per-cell Kv4 / KChIP / DPP figures clustered by type | Done | Done | Not done |
| A-current complex index and depleted-cell fractions | Done | Done | Not done (needs the per-cell depth table; `run_qc.csv` not yet transferred) |
| Per-cluster breakdown of variable subclasses | Not done (2016 types are single clusters) | Done | Not done |
| Validation against Allen RSEM counts | Not applicable | Done (ρ 0.97–0.999) | Not done (same matrix; a 5-minute step) |
| Intron/exon (nascent vs spliced) analysis | Not done; Tasic 2016 intron counts exist on GEO and would allow the same test in cortex | Done | Not done |
| Statistical tests across types (Kruskal/Mann–Whitney, BH) | Done on the preview; to be rerun on the full set | Not yet formalized beyond the regional tests | Done for the regional tests |
| Regional comparison at matched chemistry | — | Done | Done |
| Prefrontal cortex subset (PL/ILA/ORB ± ACA, 2,900–8,000 cells) | — | Not started; scope undecided | — |

Six cells of 8,137 (hippocampus) and two of 1,679 (cortex) failed quantification and are excluded; `tables/missing_cells.txt` lists them.

## Proposed experiments and analyses

In order of how much each would change the conclusions above.

1. **Kcnd2 exon/intron smFISH in intact DG and CA1 sections.** The single experiment that tests the central interpretive claim: if the spliced Kcnd2 signal is dendritic in intact granule cells while the intronic signal is nuclear in all of them, the somatic-null fraction is a localization phenotype and the “Kv4-null” language can be retired. Add a Kcnip4 or Dpp6 probe as the internal control that was normal in the same cells. Patch-seq data, which also lose dendrites, cannot answer this. Before the smFISH, a no-cost version of the same test: in laminated hippocampal RNA-seq (Farris et al. 2019, dendritic layers vs cell-body layers by laser microdissection; Cajigas et al. 2012, CA1 neuropil vs stratum pyramidale; Tushev et al. 2018 for 3′UTR isoforms), compare the exonic:intronic Kcnd2 ratio in neuropil and in tissue with the ~1:6 ratio in the SSv4 somata, normalized to Kcnip2 or Actb; a large mature-Kcnd2 excess in neuropil, with Kcnip4 and Dpp6 unchanged, supports localization. The Allen 10x whole-cell vs nuclei data from the same regions give an independent somatic-cytoplasm estimate (exon ratio whole-cell/nuclei, Kcnd2 vs Kcnip4).
2. **Validate the NP-cell complex at the protein level.** Co-labeling of Tshz2+ L5/6 cells (or Slc17a8-Cre-labeled cells, cortex only) for KChIP1 and Kv4.2 would confirm the one pyramidal class with an interneuron-type complex; it is the class predicted to be most V404M-sensitive, and it is the one where the complex composition is most unusual.
3. **V404M physiology in the two extremes within one slice.** Record A-type current inactivation in identified NP or CA1 pyramidal cells (Kv4.2-dominant) and Vip or Sncg interneurons (Kv4.3-dominant) from the same animals; the interneurons are the internal negative control. For access to NP cells, an NP enhancer AAV (Allen Institute 2025 enhancer collection) works in a knock-in background without crossing; Slc17a8-Cre covers cortical but not subicular NP cells and also labels Sncg interneurons; Rbp4-Cre and CaMKII-Cre include NP cells only within a much broader population, and NP cells have the lowest Camk2a of any pyramidal class.
4. **Confirm the KChIP-off interneuron types.** KChIP1 immunostaining or in situ in Sst HPF types (77\_Sst HPF in particular) and the Lhx6+ Lamp5 types, against the KChIP-rich Lamp5 clusters. These are candidate cells in which Kv4.2 operates without a KChIP, which changes its gating and surface expression, so they matter for V404M as well.
5. **Finish the computational set.** (a) Per-cell depth for the Allen VISp control, then the complex index and null fractions at matched chemistry for all 25 VISp subclasses. (b) Tasic 2016 intron counts from GEO to repeat the pre-mRNA-vs-spliced test in cortex. (b′) Genome alignment (STAR) of ~50 Kv4-null and ~50 Kv4-high DG cells to compare coverage across each Kcnd2 intron and exon junction: uniform intron coverage with no junction reads supports normal splicing with the mature mRNA outside the soma; a missing junction or one over-represented intron points to a splicing block. (c) Rerun the Kruskal/Mann–Whitney tests on the full Tasic set. (d) Decide the prefrontal subset (PL/ILA/ORB alone, ~2,900 cells, or with ACA, ~8,000) and run it with the pipeline as it now stands (about 6 hours per 3,000 cells on one machine).
6. **A same-chemistry check of the DPP6a ratio and Kv4.3L share.** Both are chemistry-sensitive; the Allen VISp and HPF data agree with each other, but a long-read (Iso-Seq or nanopore) library from pooled hippocampus and cortex would settle the true isoform proportions for all nine genes in one experiment and remove the chemistry caveat entirely.
7. **Low priority: 10x data.** The Allen 10x datasets are far larger but 3′-biased and cannot resolve N-terminal KChIP/DPP variants or the Kv4.3 C-terminal insert; they are useful only for the gene-level null-cell question.

## Appendix: key numbers, methods and files

**Variant nomenclature** (transcript → name): Kcnd3-202 Kv4.3S, Kcnd3-201/203 Kv4.3L; Kcnip1-203 KChIP1a, -201 KChIP1b, -207 KChIP1c, -206 KChIP1e, -202 Kcnip1-new\_short; Kcnip2-209 KChIP2a, -201/202/204 KChIP2b, -208/210 KChIP2c, -211 KChIP2x; Kcnip3-203 KChIP3a, -202 KChIP3x, -201 Kcnip3-new-cys\_rich; Kcnip4-201 KChIP4bL, -207 KChIP4b, -208 KChIP4a, -204 KChIP4e, -206 KChIP4d; Dpp6-201/206/208 DPP6a, -202/207 DPP6S, -203 DPP6K, -204/209 DPP6L; Dpp10-206/205 DPP10a, -202 DPP10c, -201/203 DPP10d. Retained-intron transcripts are dropped only when they fall in an unnamed other/fragment bucket.

**Pipeline.** Download → salmon quant (transcriptome-only vM36 index, `--validateMappings --seqBias`, `--gcBias` for paired-end) → per-transcript reads/TPM for the 9 genes → collapse. Hippocampal and VISp cells from NeMO (`.fastq.tar`, paired 50 bp); Tasic 2016 from SRA/ENA (single 100 bp). Per-cell depth = reads processed by salmon. Normalization for expression comparisons = reads per million processed (CPM). Variant fractions = pooled reads per type (figures) or per-cell fractions in cells with ≥20 reads on the family (tests). Read-weighted fractions are length-biased wherever variants have unequal annotated lengths: Kv4.3S (7 kb) vs Kv4.3L (2 kb) and the 27-kb Dpp10-202 (DPP10c) model are the affected cases; Kv4.1/4.2 shares and the KChIP and DPP6 N-terminal variants (similar lengths) change by ≤2 points between read and TPM weighting. Both datasets are dissociated whole cells (nucleus + soma; dendrites lost), not nuclei, so exonic reads report somatic mature mRNA and intronic reads report nuclear pre-mRNA.

**Pipeline lessons that affect the numbers.** NeMO allows ~3 concurrent connections per client; its tars come in two layouts (bare files or a wrapper directory); curl must be run with `-f` or HTTP error pages are saved as data. All three were fixed before the hippocampal and VISp quantifications used here; the earlier failure logs reflect those bugs, not data problems.

**Validation.** Salmon gene reads vs Allen RSEM exon counts, 8,131 hippocampal cells: Kcnd1 ρ 0.63 (Allen counts fewer), Kcnd2 0.997, Kcnd3 0.987, Kcnip1 0.975, Kcnip2 0.999, Kcnip3 0.987, Kcnip4 0.994, Dpp6 0.976, Dpp10 0.87 (salmon assigns extra reads to the 27-kb Dpp10-202 model in 2,040 cells; these are low counts). The Allen `expression_matrix.hdf5` contains intron+exon counts despite its dataset name; the exon-only CSV is the comparable one.

**Headline numbers.** Hippocampal cells by subclass: DG 2,458 · CA1-ProS 1,700 · Lamp5 891 · Vip 574 · Sst 497 · SUB-ProS 463 · Sncg 455 · CA3 322 · NP SUB 252 · Pvalb 155 · CT SUB 123. Family means over all hippocampal neurons: Kv4 86, KChIP 333, DPP 127 CPM. Somatic Kv4-null by subclass: Vip 27%, Pvalb 20%, DG 18%, CT SUB 14%, Lamp5 11%, Sst 10%, CA/SUB ≤2%. KChIP-null: Sst 19.5%, Lamp5 10.3%, others ≤1.3%. NP cells with KChIP1 >25% of KChIP: 81% (VISp), 87% (NP SUB); other pyramidal classes 0–2.5%.

**Files** (scratch directory of the analysis session; copies should be moved to the lab archive): `hip_variant_scheme_by_subclass.png`, `ctx_variant_scheme_by_type.png`, `threeway_variant_composition.png`, `hip_percell_*_genes_by_subclass.png`, `ctx_percell_*_genes_by_type.png`, `hip_acurrent_index_by_subclass.png`, `hip_lamp5_family_totals.png`, `dg_kcnd2_exon_vs_intron.png`; tables `*_variant_fracs_by_*.csv`, `*_family_zero_by_*.csv`, `hip_acurrent_index_by_subclass.csv`, `hip_exon_vs_intron_by_subclass.csv`, `threeway_variant_fracs.csv`, `threeway_tests.csv`; scripts `hip_common.py`, `variant_scheme.py`, `percell_family_genes.py`, `acurrent_index.py`, `region_compare3.py`, `run_hip_isoforms.sh`, `run_kv4_isoforms_v5.sh`, `collapse_hip.py`, `collapse_ctx.py`, `extract_markers.py`, `build_visp_runlist.py`.

**Sources.** Tasic et al. 2016, Nat Neurosci (GSE71585); Yao et al. 2021, Cell (Allen SSv4 mouse cortex + hippocampus; NeMO raw data; Allen RSEM matrices and metadata at brain-map.org); Tasic et al. 2018, Nature (VISp/ALM cells and the NP subclass; Slc17a8-IRES2-Cre); Allen Institute 2025, Cell, enhancer AAVs and transgenic lines for cortical cell types.
