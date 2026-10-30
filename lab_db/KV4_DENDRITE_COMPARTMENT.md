# Kv4 Complex mRNAs Are Somatic: Dendrite-vs-Soma Compartment Analysis

Oct 30, 2026 · @Paul

## Summary

This analysis asked whether the Kv4 channel complex mRNAs — Kcnd2 and its principal-cell partners Kcnip4, Kcnip2 and Dpp6 — are transported into dendrites in the adult mouse hippocampus, and whether any splice variant of the nine genes is preferentially dendritic. The question arose from the isoform survey (companion document, Oct 29): 18% of DG granule cells, 27% of Vip and 20% of Pvalb interneurons have no spliced Kv4 mRNA in the soma while still transcribing Kcnd2, and dendritic localization of the mRNA was offered there as the most plausible explanation.

The answer is no. In laser-microdissected hippocampal layers (Farris 2019; CA1, CA2, CA3 and DG; cell-body vs dendritic layer; 3 animals), Kcnd2, Kcnip4, Dpp6 and Kcnip2 are all depleted in the dendritic layers (log2 dendrite/soma −1.9, −1.5, −1.6 and −0.8; 12 paired comparisons, all p ≤ 0.001), while canonical dendritic mRNAs (Camk2a, Map2, Shank1, Dlg4) are enriched +1.5 to +3.2. Kcnd2 behaves like a somatic transcription-factor mRNA (Neurod6, −2.2). Moreover, the Kcnd2 signal that *is* present in the dendritic layers is 52% intronic — the same as in the cell-body layers — with intronic coverage spread uniformly over all five introns; transported mRNA would be spliced. The neuropil Kcnd2 is nuclear pre-mRNA from the interneurons and glia that reside there, not Kv4.2 mRNA in pyramidal dendrites.

The practical consequence for the survey is that **cells with low spliced Kcnd2 in the soma have low Kcnd2 mRNA, full stop**, and should be expected to have smaller A-currents unless another Kv4 gene compensates. In the survey the Kv4-null group was defined as zero reads on Kcnd1, Kcnd2 *and* Kcnd3, so for those cells (446 DG granule cells) there is no compensation at the level of detection.

A read-level follow-up on those 446 cells (STAR alignment of all Kv4-null DG granule cells plus 150 Kv4-low and 150 Kv4-high controls, unique alignments only) then asked why they lack the mRNA. Baseline Kcnd2 pre-mRNA is present in 97% of null cells and lies at the level of the non-bursting half of expressing cells (median 270 vs 1,120 intronic reads; 0.08 vs 0.14 depth per base along the gene once burst-state cells are excluded), and it is spread along the whole 500-kb gene with no pile-up in intron 1 and no hotspot at the three small annotated genes inside it. What the null cells lack is (i) the high-transcription burst state, seen in 45% of Kv4-high cells but 1% of null cells, and (ii) any mature mRNA: at matched pre-mRNA level, expressing cells carry a median 77 exonic reads and null cells 0 (98% have none). Kcnip4, Dpp6 and Kcnip2 are identical across the groups. In expressing cells the mature-mRNA level is independent of the instantaneous pre-mRNA level, so mature Kcnd2 is an accumulated, long-lived pool fed by bursts. Two models fit this snapshot — a processing/stability block specific to the null cells, or burst dynamics with non-productive baseline transcription, in which the null cells are those that have not burst within an mRNA lifetime — and separating them needs time-resolved or activity-state data (see Recommended analyses). The survey's statement that null cells have "intronic Kcnd2 identical to the other groups" rested on an Allen-derived intronic measure that does not resolve the burst state (it is uncorrelated with the cells' own exon counts, ρ = −0.03) and is corrected here. The null cells are also the low-IEG granule cells (Arc detected in 1% vs 5%, Fos 0.4% vs 4%, IEG index p < 0.001), mature but with Calb1, Prox1 and Pdzd2 reduced — a quiescent sub-state rather than a hyperactive one, which favours burst dynamics and implies that in granule cells activity promotes Kcnd2 transcription. The same Fos-low/Calb1-low profile is the ΔFosB signature of chronic hyperactivity, which matters for V404M animals (ΔFosB elevated in granule cells) and is taken up in Concerns 8.

Two splice variants do differ between compartments in a way that cell-type mixing cannot explain, and were confirmed at the level of raw splice-junction reads: the Kv4.3L cassette exon is included in 17% of Kcnd3 junction reads in the cell-body layers but 58% in the dendritic layers (p = 0.001), and the Kcnip3 cys-rich alternative acceptor falls from 19% to 4% (p = 0.002), so neuropil Kcnip3 is ~96% KChIP3a. These are the only two of ~30 testable variants that survive the controls. Their interpretation is open (see Concerns), and they are candidates, not conclusions.

The companion single-neuron dataset (Middleton 2019, 16 cultured hippocampal neurons, paired soma and dendrite) is sparse at these genes and largely GABAergic; it agrees where it can speak (Kcnip1 and Kcnd1 soma-restricted) and does not contradict anything above.

## Revision to the Oct 29 document

The statement in the companion write-up that the somatic Kv4-null cells reflect "post-transcriptional [loss], most plausibly dendritic localization of the mRNA" should be read as superseded. The data here support the alternative listed there: those cells have low cytoplasmic Kcnd2 mRNA. Whether the loss is in splicing efficiency, nuclear export or mRNA stability is not resolved by either dataset, but it is not a redistribution to dendrites.

## Datasets

| Dataset | Design | Libraries | Why it was chosen |
| --- | --- | --- | --- |
| Farris et al. 2019 (GEO GSE116342, PRJNA478217) | Adult mouse (6–8 wk) hippocampus, laser microdissection: cell-body layer vs dendritic layer (stratum radiatum for CA1/CA2/CA3; molecular layer for DG), 3 animals | 24 stranded paired-end (ISF), 42–86 M read pairs each, 70–85% salmon mapping (one CA2 dendrite library 28%) | In vivo, adult, same subfields as the SSv4 survey; deep enough for isoform quantification per layer and for genome alignment |
| Middleton et al. 2019 (GEO GSE115480, PRJNA475141) | E18 hippocampal neurons, 15 DIV, 16 neurons each micro-dissected into soma and dendrite, aRNA-amplified | 32 paired-end, 10–66 M reads, 29–81% mapping (median ~42%) | Within-neuron pairing; the authors found isoform-specific dendritic preference (alternative last exons, longer 3′UTRs) but did not look at these genes |

Four extra Farris libraries (dendritic re-runs, "-X", for CA2 and DG in animals 492 and 502) are outside the published 24-sample design and were excluded.

## Methods

1. **Isoform quantification.** salmon 1.10.3 selective alignment against the GENCODE vM36 transcriptome-only index used for the survey; variants assigned by the survey's transcript-to-isoform map. Same pipeline as the single-cell work, so variant names are directly comparable.
2. **Normalisation.** Raw TPM is unusable for compartment comparison here: 50–70% of all TPM in the Farris dendritic layers is mitochondrial (11–23% in cell-body layers), and in Middleton one unannotated rRNA-like locus (ENSMUSG00000135935) takes 20–40% of TPM in most samples, deflating every other gene in those samples by up to 3-fold. DESeq-style median-of-ratios size factors were computed within each dataset on annotated nuclear genes expressed in every sample (4,196 genes for Farris, 318 for Middleton) and applied to all TPMs. Raw and dendritic-marker-relative normalisations are retained in the tables and give the same answer for the four principal-cell genes. Variant fractions are within-gene read ratios and are unaffected by normalisation.
3. **Paired statistics.** Farris: 12 pairs (3 animals × 4 subfields), Wilcoxon signed-rank; per-subfield values (n = 3) are reported but not tested. Middleton: 16 neuron pairs.
4. **Genome alignment (Farris only).** STAR 2.7.11b on GRCm39 with the vM36 GTF; featureCounts on merged-exon and intron intervals of the 9 genes with the stranded setting matching the library (−s 1; 996 k vs 3 k exonic reads for the two strand settings confirms ISF); per-intron Kcnd2 depth by samtools bedcov on gene-region BAMs; splice-junction counts from STAR's SJ.out.tab.
5. **Cell-type composition control.** Because the dendritic layers are neuropil, they also contain glia and interneuron somata. A marker panel (Gad1/Gad2/Sst/Npy/Reln; Gfap/Aqp4/Slc1a3; Mbp/Plp1; Cx3cr1/P2ry12) was used to quantify this, and every variant shift was compared with the direction predicted by mixing pyramidal cells with interneurons using the survey's SSv4 subclass variant fractions.
6. **Junction-level validation.** For the two surviving variants, exon inclusion was recomputed from unique junction reads (inclusion junctions vs skipping junction), which is independent of the 3′UTR imbalance between transcript models that affects salmon's Kcnd3 estimates (the Kv4.3S model Kcnd3-202 alone carries a 4.8-kb 3′UTR; the Kv4.3L models end in 103 and 254 bp).

## Results, with reliability

Ratings follow the companion document: **Reliable** = within one dataset (chemistry cannot explain it), consistent across subfields or datasets, robust to the controls applied; **Tentative** = supported but resting on an interpretation the data constrain without proving, or on a single dataset with small n; **Flagged** = known or likely artifact.

### Gene level: where the complex mRNAs are

| Conclusion | Evidence | Rating |
| --- | --- | --- |
| Kcnd2 mRNA is somatic in adult hippocampal principal cells, not dendritic | log2 dendrite/soma −1.9 (normalised; −3.3 raw), 12/12 pairs down, p < 0.001, all four subfields (CA1 −2.0, CA2 −2.4, CA3 −1.3, DG −1.8); dendritic mRNAs +1.5 to +3.2 in the same samples; Neurod6 −2.2 | Reliable |
| Kcnip4, Dpp6 and Kcnip2 — the principal-cell partners — are somatic as well | −1.5, −1.6, −0.8; each p ≤ 0.001 | Reliable |
| The Kcnd2 signal in the dendritic layers is predominantly unspliced pre-mRNA | Intronic fraction of Kcnd2 reads 0.52 in dendritic vs 0.44 in cell-body layers (CA1 0.45/0.46, CA2 0.58/0.28, CA3 0.51/0.35, DG 0.54/0.66); per-intron depth per 1,000 exonic reads is flat across all five introns in both layers (45–90, no intron stands out); Kcnip4 0.76/0.82, Dpp6 0.35/0.29 | Reliable |
| The apparent dendritic enrichment of Kcnip1 (+1.5), Dpp10 (+2.3) and Kcnip3 (+0.9), and the non-depletion of Kcnd3 (+0.4), reflect interneuron (and glial) somata in the neuropil, not transport | Interneuron markers +1.9 to +3.7 and glial markers +3.0 to +4.3 in the dendritic layers; Kcnip1 and Dpp10 are 20–100× more abundant per interneuron than per pyramidal cell in SSv4 | Reliable (as an artifact assignment) |
| In cultured neurons (Middleton), Kcnip1 and Kcnd1 are soma-restricted; Kcnd2 is not significantly different | Kcnip1 −3.6 (p = 0.001), Kcnd1 −2.6 (p = 0.046), Kcnd2 −0.8 (p = 0.4); these neurons are substantially GABAergic (Gad1 670 TPM, Npy 2,680 TPM in somata) | Tentative (sparse data, 1–30 TPM for these genes) |

![Figure 1. Gene-level dendrite/soma ratios](https://github.com/ppfaff/Kv4-isoform-survey/blob/main/figures/dendrite_gene_level.png?raw=true)

**Figure 1.** log2 (dendritic layer / cell-body layer) of normalised TPM for the nine genes, per subfield (Farris) and per neuron (Middleton); points are animals or neurons. Blue Kv4, green KChIP, purple DPP.

![Figure 2. Intronic fraction by layer](https://github.com/ppfaff/Kv4-isoform-survey/blob/main/figures/dendrite_exon_intron.png?raw=true)

**Figure 2.** Intronic fraction of each gene's genome-aligned reads in cell-body (red) and dendritic (blue) layers, Farris, by subfield. Kcnd2, Kcnip4 and Dpp6 are as intronic in the neuropil as in the cell-body layer.

### Variant level: which forms move

Of the ~30 variants with enough reads to test, every compartment shift except two is either in the direction predicted by interneuron/glial admixture (KChIP2x up 3→24%; DPP6L and DPP6K up; DPP10c→DPP10d; KChIP4bL 95→70%; DPP6S 85→68%) or too thinly supported to call (KChIP1c, KChIP1e, DPP10a, KChIP4b). The two exceptions were re-tested from splice-junction reads:

| Conclusion | Evidence | Rating |
| --- | --- | --- |
| The Kv4.3L cassette exon (57 nt, chr3:105574278–105574334) is included about 3.5× more often in dendritic-layer Kcnd3 than in cell-body-layer Kcnd3 | Junction PSI 0.17 → 0.58; CA1 0.03 → 0.30, CA2 0.22 → 0.88, CA3 0.20 → 0.58, DG 0.22 → 0.29; 11/12 pairs up, p = 0.001; salmon independently gives 0.04 → 0.14 (p = 0.003). Interneurons have *less* Kv4.3L than pyramidal cells in SSv4 and glia do not express Kcnd3, so mixing predicts the opposite direction; the exon is internal coding sequence and both forms share their 3′ ends, so 3′UTR length is not involved | Tentative (robust statistically; biological reading open — see Concerns) |
| Neuropil Kcnip3 is almost entirely KChIP3a; the cys-rich alternative-acceptor form stays with the cell bodies | Junction fraction 0.19 → 0.04 (CA1 0.21 → 0.07, CA2 0.16 → 0.06, CA3 0.20 → 0.00, DG 0.25 → 0.07); 10/11 pairs down, p = 0.002. Interneurons and astrocytes carry *more* of the cys-rich form (SSv4 0.41 vs 0.28 in principal cells), so mixing predicts the opposite. The KChIP3a transcript itself is flat between layers (+0.5), the cys-rich transcript −2.9 | Tentative (robust statistically; nuclear retention is the more likely explanation) |
| Kv4.2 is a single isoform and all of it is somatic; the Kv4 pool in the neuropil shifts toward Kv4.3S/L and Kv4.1 because of interneuron admixture, not because Kv4.2 is excluded from a dendritic pool | Kv4.2 share of Kv4 reads: CA1 soma layer 0.95 → dendritic layer 0.62, with the difference made up by Kv4.3S (0.03 → 0.28), matching the interneuron Kv4 composition | Reliable |

![Figure 3. Variant composition by compartment](https://github.com/ppfaff/Kv4-isoform-survey/blob/main/figures/dendrite_variant_scheme_by_compartment.png?raw=true)

**Figure 3.** Pooled-read variant composition per dataset / subfield / compartment in the survey's colour scheme; the number beside each bar is the family read count.

![Figure 4. Junction-level test](https://github.com/ppfaff/Kv4-isoform-survey/blob/main/figures/dendrite_junction_psi.png?raw=true)

**Figure 4.** Per-animal junction-read fractions for the Kv4.3L exon (left) and the Kcnip3 cys-rich acceptor (right); each line joins the cell-body layer (left point) to the dendritic layer (right point) of one animal in one subfield.

### Read-level follow-up: Kcnd2 in Kv4-null DG granule cells

| Conclusion | Evidence | Rating |
| --- | --- | --- |
| Kv4-null granule cells have ongoing baseline Kcnd2 transcription at roughly the level of non-bursting expressing cells | Intronic reads detected (≥20) in 97.5% of 446 null cells; median 270 (IQR 156–443) vs 1,120 (328–10,687) in Kv4-high; excluding burst-state cells, mean depth per base 0.084 vs 0.139; profile flat along the gene in both | Reliable |
| The high-transcription burst state is nearly absent from null cells | >5,000 intronic reads in 45% of Kv4-high, 10% of Kv4-low, 1.1% of Kv4-null cells | Reliable |
| Null cells have no mature Kcnd2 mRNA even at matched pre-mRNA level | At <450 intronic reads: Kv4-high n = 49, median 77 exonic reads, 0% zero; Kv4-null n = 337, median 0, 98% zero; 93% of null cells have no splice-junction read | Reliable |
| Mature Kcnd2 mRNA in expressing cells is an accumulated pool, not a readout of current transcription | Exonic reads flat (72–125) across intronic bins from <150 to >5,000 | Reliable |
| The Kcnd2 deficit is gene-specific | Kcnip4 (26k/23k/23k intronic), Dpp6 (12k/10k/10k) and Kcnip2 (286/219/245) identical across groups; exon counts likewise | Reliable |
| The intronic signal is Kcnd2 pre-mRNA, not an intronic lncRNA/TEC | 10-kb windows: no hotspot at Gm42583 (21.26 Mb), ENSMUSG00000122200 (21.43 Mb, antisense) or Gm42584; depth uniform across intron 1 and introns 2–5 | Reliable (SSv4 is unstranded, so antisense cannot be excluded by strand, only by position) |
| Mechanism: processing block vs burst dynamics | Not separable from a snapshot (see Concerns 7) | Open |
| 14 "null" cells (3%) with >1,500 intronic reads do carry 30–50 exonic reads by STAR | Salmon/STAR discordance at the zero boundary; these are expressing cells misclassified by salmon's zero | Flagged (minor) |

![Figure 5. Read-level Kcnd2 in DG granule cells](https://github.com/ppfaff/Kv4-isoform-survey/blob/main/figures/dg_kcnd2_readlevel.png?raw=true)

**Figure 5.** A: per cell, Kcnd2 intronic vs exonic unique reads, coloured by the survey's Kv4 group; the dashed line marks the burst state. B: median exonic reads per intronic-read bin. C: mean depth per base per cell in 10-kb windows along the gene, burst-state cells excluded (exons shaded).

### Activity state of the Kv4-null granule cells

Kv4.2 is activity-downregulated (Kcnd2 mRNA falls in DG after kainate; seizures load miR-324-5p onto Ago2 and drive Kv4.2 internalisation), so if the null state were an activity state the expectation was high-IEG cells. The Allen exon counts for all 2,458 DG cells say the opposite.

| Conclusion | Evidence | Rating |
| --- | --- | --- |
| Kv4-null granule cells are the low-IEG (quiet) cells | Detection of Arc 1.1% vs 5.0% (Kv4-high), Fos 0.4% vs 4.3%, Egr1 17% vs 24%, Egr4 16% vs 29%, Dusp1 0.7% vs 5.3%, Homer1 48% vs 59%, Fosl2 16% vs 25%; 16-gene IEG index median 1.07 vs 1.22 (all p ≤ 0.001, Mann–Whitney); Bdnf, Kcnip4, C1ql2 unchanged; sequencing depth and library complexity identical; the difference holds within every library-size tertile | Reliable |
| Null cells are mature but express several mature-granule-cell transcripts at lower levels | Dcx, Calb2, Sox11, Sox4, Igfbpl1 identical (not adult-born); Calb1 −24%, Prox1 −29%, Pdzd2 −27%, Dsp and Tubb3 lower (p ≤ 0.001); mild shift toward Allen clusters 363/364_DG | Reliable (within dataset) |
| Interpretation: a quiescent sub-state of mature granule cells in which Kcnd2 transcription bursts are rare and the mRNA pool has decayed | Consistent with sparse granule-cell firing; implies activity *promotes* Kcnd2 bursting in granule cells rather than suppressing it, opposite to the seizure literature at the protein level | Tentative (correlational snapshot; direction of causality untested) |
| The Fos-low / Calb1-low combination is also the ΔFosB-high signature of chronically hyperactive granule cells (Chin lab) | Not the parsimonious reading in wild-type Allen mice (Fosb mRNA itself is lowest in null cells, 0.7% vs 3%), but directly relevant to V404M animals, where ΔFosB is elevated in granule cells: a V404M dentate could carry Kv4-null cells of either kind, distinguishable only by ΔFosB itself | Open |

ΔFosB could not be read from the SSv4 junctions (Fosb is detected in ~3% of cells; no Fosb-locus junction reads in the 746 aligned cells).

![Figure 6. Activity-regulated genes by Kv4 group](https://github.com/ppfaff/Kv4-isoform-survey/blob/main/figures/dg_activity.png?raw=true)

**Figure 6.** Left: mean log2(CPM+1) of activity-regulated and granule-cell genes in the three Kv4 groups (Allen exon counts, 2,458 DG cells). Right: per-cell IEG index.

### Technical findings worth recording

| Finding | Detail | Rating |
| --- | --- | --- |
| Salmon underestimates Kv4.3L in bulk data by ~4-fold | Salmon cell-body-layer Kv4.3L fraction 0.04 vs junction PSI 0.17; the 4.8-kb 3′UTR on the Kv4.3S model alone pulls shared reads to S. This is the same transcript-model asymmetry flagged in the companion document for SSv4, where the Kv4.3S/3L balance was already the one unsettled Kv4 number. The junction PSI is the number to trust; the survey's Kv4.3L fractions should be regarded as lower bounds | Flagged (survey value), Reliable (junction value) |
| DPP6a is 0.2% of Dpp6 in the Farris cell-body layers but ~35% of Dpp6 in SSv4 pyramidal cells | Not pursued; most likely a 5′-end coverage difference between the two chemistries acting on alternative-first-exon transcripts. Reinforces the rule that isoform fractions are only compared within a dataset | Flagged |
| Raw TPM compartment ratios are dominated by mitochondrial and rRNA composition | See Methods 2; all marker genes, including Camk2a and Map2, appeared dendritically depleted before correction | Flagged (method) |
| The Farris libraries are forward-stranded (ISF), not the dUTP (ISR) layout the kit name suggests; ~half the Middleton libraries are called unstranded by salmon | Harmless for quantification (salmon infers per library) but required for intron counting; recorded so the strand setting is not re-derived wrongly | Reliable |

## Concerns and open questions

1. **Laminar interneurons.** The neuropil layers hold specific interneuron populations (CCK/Sncg basket and Schaffer-associated cells, neurogliaform Lamp5 cells, O-LM and other Sst cells in stratum oriens/lacunosum-moleculare), and SSv4 classifies cells by transcriptome, not position. The mixing control used subclass averages; if the neuropil-resident subset of an interneuron subclass has a different Kv4.3L (or Kcnip3) splice ratio from the subclass average, the Kv4.3L result could still be a cell-identity effect rather than a localisation effect. Nothing in the present data excludes this.
2. **Kv4.3L: localisation vs identity.** If the exon-inclusion difference is genuinely compartmental, it implies either preferential transport of the L-form mRNA or compartment-specific splicing — both unusual for an internal cassette exon, and the 19-residue L insert (C-terminal, carrying a PKC consensus site) would be a plausible thing for a dendrite to control. But the simpler reading remains that the cells whose Kcnd3 fills the neuropil are L-rich. Kcnd3 at the gene level is not depleted in the neuropil (+0.4), so a substantial Kcnd3 mRNA pool is there in either case.
3. **Kcnip3 cys-rich form: retention vs export.** Its share declines monotonically from nucleus-rich SSv4 cells (28%) through cell-body layers (19%) to neuropil (4%), the pattern of a nuclear-retained or inefficiently spliced species produced by a frame-shifting 4-nt alternative acceptor. Whether it is ever exported and translated is unknown; the KChIP3a form is the only one with a plausible dendritic presence. Given that KChIP3 is also DREAM, a nuclear Ca²⁺-regulated transcriptional repressor, the N-terminal variant question for this gene is at least as much about protein targeting as about mRNA, and the mRNA data cannot speak to that.
4. **Sensitivity floor.** Bulk laminar data cannot exclude a small dendritic Kcnd2 pool (tens of copies per dendrite) of the kind reported in cultured neurons (Lee et al. 2011; Gross et al. 2011). The conclusion is that no *substantial* fraction of Kcnd2 mRNA is dendritic in vivo. The direct FMRP–Kv4 protein interaction (Zhan et al. 2024) explains the Fmr1 A-current phenotype without a dendritic mRNA pool, which is consistent with these data.
5. **Only one in vivo dataset.** Every Reliable rating above rests on Farris 2019 (3 animals). The same-direction result in all four subfields is strong internal replication, but an independent preparation has not been examined.
6. **Middleton is weak at these genes.** The 32 single-neuron libraries are sparse (318 genes above 5 TPM in all samples), the neurons are largely GABAergic, and half are called unstranded. It should not be used for the auxiliary-subunit questions beyond the two significant somatic calls.
7. **Why the null cells lack mRNA is not settled, but the activity data favour burst dynamics.** A processing/stability block and a burst-dynamics model (long-lived mRNA fed by bursts; non-productive baseline transcription in all cells; null cells = cells that have not burst within an mRNA lifetime) give the same read-level snapshot. The null cells being the low-IEG population, with several mature-granule transcripts also reduced, fits a quiescent state in which bursts have not occurred recently better than a cell-specific processing defect, but this is correlational. Metabolic labelling (SLAM-seq / TimeLapse-seq on DG) or smFISH with intron and exon probes in the same section would settle it. A miRNA acting on the accumulated pool (miR-324-5p is translational, not decay, per Gross 2016) is not excluded but is not needed.
8. **ΔFosB and V404M.** In V404M animals ΔFosB is elevated in granule cells (lab observation), the marker of chronic dentate hyperactivity, and ΔFosB represses c-Fos and Calb1. A Kv4-null granule cell in a V404M dentate could therefore be either a quiet cell (as in wild type here) or a chronically overactive, ΔFosB-high cell with the same Fos-low / Calb1-low mRNA profile. Whether V404M changes the size or nature of the Kv4-null population is a question the existing data cannot address.

## Recommended additional analyses

1. **Compensation check in Kv4-low cells (immediate, from existing tables).** Among SSv4 principal cells with *low but non-zero* Kcnd2, test whether Kcnd1 or Kcnd3 is elevated relative to Kcnd2-high cells of the same type. The null group (zero on all three Kcnd genes) already shows none; the low group is the one that matters for predicting graded A-current loss. This is a few lines against `hip_isoform_per_cell_annotated.csv` on Lambda.
2. **Independent in vivo replication (rat).** Tushev et al. 2018 (3′-end sequencing, CA1 neuropil vs somata) and Perez et al. 2021 (neuropil subcellular fractions) provide a second species and a different library design. Re-running the Kcnd2 exon/intron and Kv4.3L junction tests there would need a rat STAR index (~2 h on Lambda) and a one-time nomenclature mapping of the rat transcripts. This is the single most useful next step for the two tentative variant results.
3. **Interneuron junction-level Kv4.3L (SSv4).** Align a few hundred hippocampal interneurons from the Yao 2021 SSv4 set (Lamp5, Sncg, Sst, Pvalb, Vip; ~50 per subclass) with STAR and count the Kv4.3L inclusion junctions per subclass. This replaces the salmon-based L fractions (which this work shows are 4-fold low) with a direct measurement, and tests whether any interneuron subclass is L-rich enough to account for the neuropil PSI of 0.58–0.88.
4. **Kcnip3 cys-rich acceptor: nuclear vs cytoplasmic.** Any nuclear/cytoplasmic fractionation RNA-seq of mouse brain (several exist) would show directly whether the 4-nt acceptor form is nuclear-retained. If it is, the cys-rich variant can be dropped from the functional variant list.
5. **smFISH for Kcnd2 exon vs intron in DG (unchanged from the companion document).** The present result sharpens its prediction: Kcnd2-exon puncta should be somatic and should be reduced or absent in a ~18% subset of mature granule cells, while intronic puncta remain nuclear in all cells.
6. **Kv4.3L in CA2.** CA2 stratum radiatum Kcnd3 is 88% L-form. If the Kv4.3L result survives (2) and (3), CA2 is where to look first, and the CA2 neuropil interneuron population is well characterised.
7. **Burst dynamics vs processing block in the Kv4-null granule cells.** The IEG comparison (done; null cells are IEG-low) favours burst dynamics. Decisive data: SLAM-seq or TimeLapse-seq on dentate gyrus to measure Kcnd2 mRNA half-life and the productive fraction of baseline transcription; dual intron/exon smFISH (item 5) distinguishes the two if nuclear intronic foci are present in exon-negative cells. A direct test of the activity link: Kcnd2 exon smFISH with Fos after a controlled activity manipulation (e.g. novel-environment exposure vs home cage), asking whether the Kv4-null fraction shrinks in the active dentate.
8. **V404M dentate.** smFISH (or snRNA-seq) on V404M vs wild-type dentate for Kcnd2 exon/intron together with Fos, Calb1 and ΔFosB (antibody or Fosb-ΔFosB-specific probe), to ask whether the Kv4-null fraction changes and whether V404M null cells are quiet or ΔFosB-high.

## Implications for Kv4.2-V404M

The companion document ranked cell types by exposure to the mutation using somatic Kv4.2 share and complex completeness. This analysis removes the main caveat on that ranking: the 18% of DG granule cells (and the Vip and Pvalb subsets) without somatic Kv4 mRNA are not cells that have shipped their Kv4.2 to the dendrites; they are cells with little Kv4 mRNA, and within the null group no Kv4.1 or Kv4.3 compensation is seen. Such cells should have smaller A-currents, and in a V404M heterozygote they also contribute less mutant channel — so the mutation's effect on a circuit is carried by the Kv4-expressing majority, and the Kv4-low minority is a pre-existing source of A-current heterogeneity that the mutation does not create. The one place a dendritic variant question remains live is Kv4.3, which is not the V404M gene.

## Data, code and provenance

- Code: `pipeline/build_dendrite_runlist.py`, `run_dendrite_isoforms.sh`, `run_dendrite_star.sh`, `make_intron_saf.py`, `collapse_dendrite.py`; `scripts/dendrite_compartment.py`, `scripts/dendrite_junctions.py`; read-level DG follow-up `pipeline/build_dg_runlist.py`, `run_dg_star.sh`, `dg_kcnd2_windows.sh`, `scripts/dg_readlevel.py`, `scripts/dg_activity.py` (activity genes via `pipeline/extract_markers.py` with `GENES=`) in https://github.com/ppfaff/Kv4-isoform-survey.
- Tables: `tables/dg_activity_by_group.csv`, `tables/dg_kcnd2_readlevel_summary.csv`, `tables/dg_kcnd2_per_cell_q255.csv`, `tables/kcnd2_locus_annotation.tsv` (DG read-level); `tables/dendrite_*.csv` (gene level, variant fractions, marker QC, exon/intron, junction PSI, variant-vs-mixing), `tables/star_*.tsv`, `tables/sj_kcnd3_kcnip3.tsv`.
- Raw quantifications: `~/kv4_dendrite/` on PFAFF-LAMBDA (salmon quant for 56 libraries; STAR region BAMs, SJ.out.tab and featureCounts output for the 24 Farris libraries).
- Reference: GENCODE vM36 / GRCm39 throughout; variant nomenclature as in the companion document.
- Statistics are as of 2026-10-30.
