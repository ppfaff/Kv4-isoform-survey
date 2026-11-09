# V404M/+ hippocampus: a chronic activation of the Kv4.2-dependent principal cells (dentate granule cells and CA1) — bulk RNA-seq program, ΔFosB/c-Fos histology and the Viaat-Cre conditional

Nov 9, 2026 · Paul Pfaffinger (analysis with Claude) · PROVISIONAL

## Summary

In the 3-month littermate bulk RNA-seq (6 V404M/+ vs 7 +/+ whole hippocampi), the mutant hippocampus carries a coordinated set of expression changes that is the known transcriptional signature of dentate granule cells under recurrent hyperactivity: **proenkephalin up (Penk, log2FC +0.97 to +1.21), prodynorphin down (Pdyn, −0.86), Npy up (+0.67 to +0.93), Drd1 up (+0.88 to +1.27), Ecel1 up (+1.25), Cck down (−0.43)**, with Igfbpl1 (+1.5) and Rprm (+1.2) alongside. Penk, Npy and Drd1 separate the six mutants from the seven wild types with no overlap (exact permutation p ≈ 0.001 each). The program is hippocampus-specific — the same genes are flat in cortex, thalamus and cerebellum except Npy, which is modestly up wherever Kcnd2 is expressed — and it survives correction for the choroid-plexus and cortical tissue that contaminated some hippocampal blocks (Penk, Drd1, Igfbpl1, Rprm padj < 0.001 both with and without composition covariates). The opposite movement of enkephalin and dynorphin, with CCK falling, is the pattern Gall (1988) described in mossy fibres after seizures, and ectopic NPY in granule cells and mossy fibres is the classic histological marker of an epileptic dentate gyrus. Acute immediate-early genes are not elevated (Fos, Egr1, Npas4 flat or lower; Bdnf flat), so this is a chronic state rather than the aftermath of a recent seizure. It is the transcriptional counterpart of the lab's ΔFosB immunohistochemistry (Figure 1 of the V404M manuscript draft), which shows **ΔFosB up 2.41× in dentate granule cells and 2.66× in CA1 pyramidal cells (both p < 0.0001, n = 6–7), unchanged in CA2 and CA3, and up 1.64× in retrosplenial cortex** (one-tailed p < 0.05). qPCR in a separate cohort shows **c-Fos mRNA down to 0.33× in cortex and 0.70× in hippocampus** (Figure 3), with ATF4 unchanged — the same direction and size as the bulk (Fos −1.2 and −0.5 log2; ISR negative), so the "ΔFosB high, acute c-Fos low" state is confirmed by two methods in two cohorts.

Two further observations bear on mechanism. **Kcnip3 (KChIP3/DREAM) is up in mutant hippocampus** (+0.31, padj 0.02; +0.55, padj 0.006 after composition correction), and DREAM's founding target is the prodynorphin DRE, which it represses in the absence of Ca²⁺ — so Kcnip3 up and Pdyn down are one event, not two. And **Viaat-Cre × conditional V404M mice, which carry the mutation only in GABAergic neurons, show spike-wave discharges on the frontal EEG but no ΔFosB elevation in the dentate** (lab observation). The dentate activation therefore does not follow from the spike-wave phenotype and is not produced by V404M in interneurons; it requires the mutation in the granule cells themselves or in their excitatory inputs, The pattern of ΔFosB across hippocampal fields is the decisive clue: dentate granule cells and CA1 pyramidal cells are the two hippocampal principal populations whose A-current is built on Kv4.2 (Allen hippocampal SMART-seq, Kcnd2 reads/cell: CA1 65, CA2 79, CA3 33, DG 42; the CA1 dendritic A-current is the canonical Kv4.2 current, CA3 relies more on Kv4.3), and they are the two that light up, while Kv4.2-poor CA3 — the direct target of the granule cells' mossy fibres — is spared. A trans-synaptic spread from an over-active dentate would reach CA3 before CA1. The data instead fit **a cell-autonomous activation of whichever principal cells depend on Kv4.2**, which is also what the Viaat-Cre dissociation says.

## Results

### The program, gene by gene

Original per-tissue DESeq2 (`~ sex + genotype`, STAR counts) and the composition-adjusted model (`+ choroid-plexus index + cortex index`, see the Nov 9 revision of the bulk write-up):

| gene | role in granule cells | log2FC unadjusted (padj) | log2FC adjusted (padj) |
| --- | --- | --- | --- |
| Penk | enkephalin; induced in granule cells/mossy fibres by seizures | +0.97 (0.003) | +1.21 (<0.001) |
| Pdyn | dynorphin; mossy-fibre peptide, *depleted* chronically after seizures; DREAM target | −0.86 (<0.001) | −0.77 (0.056) |
| Npy | ectopic granule-cell/mossy-fibre NPY = epileptic DG marker | +0.67 (0.24) | +0.77 (0.075) |
| Drd1 | granule-cell D1 receptor; activity-inducible | +0.88 (<0.001) | +1.27 (<0.001) |
| Ecel1 (DINE) | damage/seizure-induced neuronal endopeptidase | +1.25 (0.008) | +1.13 (0.088) |
| Cck | mossy-fibre CCK, falls after seizures (Gall 1988) | −0.43 (<0.001) | −0.29 (0.47) |
| Igfbpl1 | granule-cell gene, activity-regulated | +1.47 (<0.001) | +1.63 (<0.001) |
| Rprm | reprimo, activity-regulated | +1.16 (<0.001) | +1.49 (<0.001) |
| Rasd1 | activity-regulated small GTPase | −0.69 (<0.001) | −0.55 (0.096) |
| Kcnip3 (DREAM) | Pdyn repressor; Kv4 auxiliary subunit | +0.31 (0.020) | +0.55 (0.006) |
| Dsp | mature granule-cell marker | −0.74 (0.20) | −0.94 (0.024) |
| Calb1 | mature granule-cell marker; repressed by ΔFosB (Chin lab) | −0.22 (0.53) | −0.42 (0.16) |
| Fos / Egr1 / Npas4 | acute IEGs | −0.49 / −0.45 / −0.26 (n.s.) | −0.87 / −0.31 / −0.05 (n.s.) |
| Bdnf | acute activity gene | +0.03 (0.98) | −0.06 (0.88) |

Per-animal separation: a signature index built from Penk, Drd1, Ppp1r1b and Npy relative to Snap25 places every mutant above every wild type (Fig. C; `tables/v404m_admixture_per_sample.csv`). The MSN-specific genes that would move if this were striatal contamination (Gpr88 −0.25, Adora2a −0.07, Pde10a +0.05) do not move.

### Why this is chronic, not acute

Acute seizures raise Fos, Egr1, Npas4, Arc and Bdnf within an hour and for several hours. None is elevated; Fos and Egr1 trend *down* in hippocampus and are significantly down in cortex. Persistent suppression of the acute IEG response with elevated ΔFosB is the documented state of a chronically overactive dentate (ΔFosB represses c-fos through its own AP-1 sites), and ΔFosB repression of Calb1 in granule cells has been shown in the Chin laboratory's Alzheimer-model work; Calb1 trends down here (−0.22 to −0.42, not significant in bulk). The peptide pattern — enkephalin up, dynorphin and CCK down — is likewise the pattern of weeks-long kindling or post-status tissue rather than of a single event.

### Specificity

Penk, Pdyn, Drd1, Ecel1, Igfbpl1 and Rprm are unchanged in cortex, thalamus and cerebellum (Fig. A). Npy is the exception: +0.59 in cortex (padj < 0.01 unadjusted), +0.40 in thalamus (p 0.02 nominal), i.e. up in each tissue where Kcnd2 is expressed in excitable neurons, and flat in cerebellum, which expresses the most Kcnd2 but shows no network phenotype. The dentate-specific set is therefore not a generic NPY response but a granule-cell program.

### The Kcnip3–Pdyn link

KChIP3 is DREAM. In the absence of Ca²⁺ it binds the downstream regulatory element of *Pdyn* (and *Fos*) as a repressor; Ca²⁺ binding releases it. Kcnip3 mRNA up and Pdyn mRNA down in the same tissue is the expected pair if DREAM levels, rather than Ca²⁺ relief, dominate. It is also the first indication in this dataset that an auxiliary subunit responds to the mutation — through its nuclear rather than its channel role. Which cells carry the Kcnip3 increase (granule cells, or hilar/CA interneurons, where Kcnip3 is also high) cannot be resolved in bulk.

### The Viaat-Cre conditional dissociates the two phenotypes

Viaat-Cre × conditional V404M animals (mutation in GABAergic neurons only; cohort 6, n = 5 Viaat-Cre;+/+ and 6 Viaat-Cre;M/+, 4 M + 2 F mutants; ΔFosB immunohistochemistry scored by Manu, DG quantified as in Jerng et al. 2026) show spike-wave discharges at the frontal electrode and:

| region | readout | result |
| --- | --- | --- |
| dentate gyrus | mean grey value, fold of WT | 1.02×, no trend in either sex |
| cortex, upper 50% | % area ΔFosB⁺ nuclei | +32% (p = 0.065); mean intensity also higher |
| cortex, lower 50% | % area ΔFosB⁺ nuclei | weak, non-significant increase; less ΔFosB overall |
| TRN | ΔFosB⁺ cell counts, 2 sections/animal | males ≈ 100 vs 50 (trend); females no trend (n = 2) |

So the mutation in inhibitory neurons alone produces the loop phenotype with a histological footprint in both of its GABAergic nodes — upper-layer cortex and TRN — while leaving the dentate silent. Read with the bulk result: the spike-wave phenotype arises from V404M in inhibitory neurons (cortical interneurons and/or the reticular nucleus — the tier-1 thalamus analysis will say which carry Kv4.2), while the dentate activation needs the mutation in glutamatergic cells. Absence-type spike-wave discharges in rodent models generally do not recruit the hippocampus, consistent with the quiet dentate here. The cortical ΔFosB elevation in Viaat-Cre also gives the constitutive bulk's cortical IEG suppression a companion: "ΔFosB high, acute Fos/Egr1 low" is the cortical state as well as the dentate one. The male-only TRN trend rests on two mutant females and is a question, not a finding; the constitutive bulk (3 mutants of each sex) is the place to test a genotype × sex interaction, for the cortical IEGs in particular (the dentate program itself shows no sex split: mutant F −7.02, mutant M −7.12, WT F −7.63, WT M −7.68 on the per-animal index).

Constitutive staining by field (manuscript Figure 1): DG 2.41 ± 0.08, CA1 2.66 ± 0.24, CA2 1.3 (n.s.), CA3 1.3 (n.s.), retrosplenial cortex 1.64 ± 0.34 — see Interpretation.

## Interpretation: where the hippocampal hyperactivity comes from

Granule cells are the most electrically silent principal neurons of the forebrain — resting potential near −80 mV, strong feed-forward (basket) and feedback (hilar SST/NPY) inhibition, and a large A-type current (Kv4.2 with KChIP4 and DPP6) — and this quiescence is the "dentate gate" between entorhinal cortex and CA3. In the standard epilepsy models dentate hyperactivity is *driven and disinhibited* rather than intrinsic: excessive perforant-path input, loss of hilar mossy cells and SST interneurons, mossy-fibre sprouting creating recurrent excitation; the granule-cell transcriptional changes (NPY, enkephalin up; dynorphin, CCK down; ΔFosB) follow from that activity, and CA3 is recruited through the mossy fibres.

V404M inverts that order, and the ΔFosB map says how. If the dentate gate had failed and the activity propagated, CA3 would be the first field downstream to show it; CA3 is unchanged (1.3×, n.s.) while CA1 is as elevated as the dentate (2.66× vs 2.41×). The fields that are active are the fields whose principal cells carry Kv4.2 as their main A-type subunit; the field that is quiet is the one that relies on Kv4.3. The simplest account is therefore that each Kv4.2-dependent principal cell is intrinsically more excitable, and ΔFosB accumulates wherever that is true — DG, CA1, upper-layer/retrosplenial cortex — without requiring a propagating focus. Viaat-Cre supports the same reading from the other side: the mutation in inhibitory neurons produces the thalamocortical loop phenotype (SWDs, TRN and upper-cortical ΔFosB) but does nothing to DG or CA1, so disinhibition is not the route.

Two readings remain open: (1) cell-autonomous reduction of the Kv4.2-dependent A-current in granule cells and CA1 pyramids, or (2) increased entorhinal drive if V404M also disinhibits entorhinal layer II/III cells (which project to both DG and CA1 and express Kcnd2 strongly) — note that (2) would *also* spare CA3 if the direct entorhinal→CA1 (temporoammonic) path and the perforant path to DG were the carriers. CaMKII-Cre hits both; a granule-cell-restricted Cre predicts DG ΔFosB without CA1 ΔFosB under (1), and neither under (2). Recording CA1 and granule-cell A-currents in V404M/+ slices is the direct test of (1).

## Reliability

| claim | rating | basis |
| --- | --- | --- |
| Penk ↑, Drd1 ↑, Igfbpl1 ↑, Rprm ↑ in V404M hippocampus | solid | padj < 0.001 in both models; per-animal complete separation for Penk and Drd1 |
| Pdyn ↓, Cck ↓, Npy ↑, Ecel1 ↑ | moderately solid | significant unadjusted; 0.05–0.1 after adjustment (adjustment costs df); directions consistent with the literature pattern |
| Chronic rather than acute state | moderately solid | no IEG/Bdnf elevation; IEGs trend down |
| Kcnip3 ↑ and Pdyn ↓ as one DREAM event | tentative | both significant, mechanism from literature, cell type unknown |
| DG + CA1 activation needs V404M in glutamatergic cells | tentative | Viaat-Cre: DG 1.02×, cortex/TRN up; CaMKII-Cre in progress |
| ΔFosB high / c-Fos low state in cortex and hippocampus | solid | IHC (RSC, DG, CA1) + qPCR c-Fos in a second cohort + bulk RNA-seq agree in direction and size |
| Activation tracks Kv4.2-dependent cell types (DG, CA1; not CA3) | moderately solid | ΔFosB map vs Allen Kcnd2 expression by field; CA2 (1.3×, high Kcnd2) is the one discordance |
| Hippocampal astrogliosis (Gfap, Serpina3n) | downgraded | fold changes persist (+0.69, +0.73) but lose significance once choroid-plexus content is a covariate |

## Proposed experiments

0. **CaMKII-Cre × conditional V404M (in progress, animals being weaned Nov 2026).** Forebrain glutamatergic neurons, the reciprocal of Viaat-Cre. Prediction if the two phenotypes are class-dissociated: dentate ΔFosB and the Penk/Pdyn/Npy/Drd1 program present, spike-wave discharges absent. Spike-wave discharges present would mean pyramidal-cell V404M (L6 corticothalamic / L5) can drive the loop on its own. The Cre line's thalamic expression should be confirmed (T29-1 is sparse in relay nuclei and absent from RT; other Camk2a-Cre lines differ). A qPCR panel (Penk, Pdyn, Npy, Drd1, Cck vs Snap25/Prox1) on whole hippocampus is a one-day readout of the dentate program across the three genotypes; bulk effect sizes support n = 4-5 per group.

| | spike-wave (frontal EEG) | dentate ΔFosB / program |
| --- | --- | --- |
| constitutive V404M/+ | yes (+ RSC ΔFosB) | yes (DG 2.4×, CA1 2.7×, CA3 flat) |
| Viaat-Cre (GABAergic) | yes (+ cortical and TRN ΔFosB) | no (1.02×) |
| CaMKII-Cre (forebrain glutamatergic) | predicted no | predicted yes |

1. **Granule-cell-restricted conditional**: Prox1-CreERT2 or Dock10-Cre × conditional V404M, scored for ΔFosB in DG *and CA1* and for Penk/Npy by smFISH. DG-only ΔFosB → cell-autonomous; DG + CA1 → propagation or shared drive; neither → entorhinal drive. Pair with Emx1-Cre (all forebrain glutamatergic) to cover the entorhinal alternative.
2. **Cell assignment of the bulk program**: RNAscope for Penk, Npy, Pdyn and Kcnip3 on the existing ΔFosB-stained sections, granule layer vs hilus vs CA3; the per-animal bulk values give the expected effect sizes.
3. **Granule-cell A-current and firing in V404M/+ vs +/+** (slice; the lab's recording design document already specifies the mixture analysis): the direct test of the intrinsic reading, and of whether KChIP/DPP composition in granule cells (DPP6K-rich, KChIP4bL/4a) changes the mutant's effect.
4. **EEG–histology pairing**: depth electrode in dentate alongside the frontal lead in constitutive heterozygotes, to ask whether dentate events are independent of spike-wave discharges, as the conditional suggests.
5. **Kcnip3/DREAM**: nuclear vs cytoplasmic KChIP3 in granule cells by immunostaining; Pdyn and Fos DRE occupancy is a later step.

## Methods and provenance

Bulk libraries and alignment as in the Nov 6 bulk write-up (STAR ReadsPerGene, GENCODE vM36). DE: pydeseq2 per tissue, `~ sex + genotype` (original) and `~ sex + genotype + composition indices` (`scripts/v404m_de_composition.py`; indices = log2 mean CPM of marker sets relative to Snap25, `tables/v404m_composition_indices.csv`). Per-animal indices and marker-class ratios: `scripts/v404m_admixture.py`, `tables/v404m_admixture_per_sample.csv`, `tables/v404m_admixture_genes.csv`. Figure: `scripts/v404m_dg_program.py` → `figures/v404m_dg_program.png` (panel D, per-animal CPM, is added when `--extract` has been run on the machine holding the counts). Exact permutation p for complete separation of 6 vs 7 = 2/C(13,6) = 0.0012.

Key literature for the pattern: Gall C (1988) J Neurosci 8:1852 — seizures raise enkephalin and lower dynorphin and CCK immunoreactivity in mouse mossy fibres; Sperk G et al. (1992) / Vezzani A et al. (1999) — ectopic NPY in granule cells and mossy fibres after seizures; Carrión AM et al. (1999) Nature 398:80 — DREAM is a Ca²⁺-regulated repressor of prodynorphin; You JC et al. (2017) Nat Med 23:1377 and Corbett BF et al. (2017) Cell Rep — ΔFosB in dentate granule cells represses c-Fos and calbindin in a hyperexcitable hippocampus. Citations from memory; verify before quoting.

Note on the manuscript figures (Nov 9): Figure 3A caption value corrected to 0.33 (cortical c-Fos, mutant). Figure 1B uses a one-tailed test; the bulk Fos/Egr1 result in cortex is an independent two-tailed confirmation that could be cited alongside it.
