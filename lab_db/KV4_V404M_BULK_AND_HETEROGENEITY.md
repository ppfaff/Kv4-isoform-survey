# Kv4.2-V404M Bulk RNA-seq, and Within-Type Heterogeneity of Kv4-Complex Composition

Nov 6, 2026 · @Paul

## Summary

Two separate pieces of work, reported together because the second changes how the first should be read.

**V404M bulk RNA-seq** (48 libraries, littermate pairs, 4 tissues, 3-month-old mice) answers the four questions put to it, three of them negatively. Mutant animals are heterozygous with **both alleles expressed at parity** (knock-in fraction 0.495; the two haplotype windows carry 5,510 and 5,485 reads), so the missense change does not destabilise its own transcript. There is **no integrated stress response, no unfolded-protein response and no apoptotic program** in any tissue: the XBP1/IRE1 splice shows no mutant elevation anywhere, and the one nominal hit (thalamus, p = 0.022) runs the wrong way. The **Kcnd2 exonic:intronic ratio is unchanged**, and the direction in hippocampus is *less* intron retention in mutants (0.0292 vs 0.0390, p = 0.0221 nominal, q ≈ 0.8 across ~36 gene×tissue tests) — the opposite of an increased unprocessed/null fraction. The **ΔFosB junction** is elevated 2.3-fold in mutant hippocampus (0.0667 vs 0.0294, p = 0.100), the largest of the four tissues, with cortex flat and cerebellum at zero; the tissue-specificity matches the lab's ΔFosB immunostaining, but at 109 reads pooled across 48 libraries this is underpowered, and ΔFosB accumulates through protein half-life rather than transcript abundance, so bulk junction counts are the wrong assay for it either way.

What the mutants *do* show is **a reactive astrocyte state in hippocampus and cortex and nowhere else** (Gfap +0.79 and +0.62; Serpina3n +0.89 and +0.54; Timp1, Ptx3, Vim, Osmr in hippocampus). Three findings constrain its interpretation. It is **not proliferation**: astrocyte identity genes (Aqp4, Slc1a3, Gja1, Aldoc) are flat at −0.04 in hippocampus. It is **not the complement/neurotoxic program**: C3 is *down*. And **per-astrocyte homeostatic capacity is preserved** — Kcnj10 (Kir4.1) +0.17 in hippocampus, Slc1a2 (GLT-1) −0.01, so the two canonical routes by which reactive astrocytes amplify hyperexcitability are not engaged. Decisively, the regional pattern does **not** track mutant protein dose: cerebellum expresses the most Kcnd2 of any tissue (341 CPM vs hippocampus 166, thalamus 112, cortex 101) and shows no reactive signature at all. It tracks which circuits are epileptogenic. On this evidence the astrogliosis is a **readout of circuit dysfunction, not a driver of it**.

**Within-type heterogeneity.** Asked whether the large per-cell variance in Kv4/KChIP/DPP composition reflects lumping distinct neuron types or intrinsic variation within coherent ones, the answer is **intrinsic**. Every subclass in all three single-cell datasets shows within-subclass spread 4–25× the per-cell measurement noise floor, and the ratios scale as √sequencing-depth across the three (observed 2.02, 3.13, 1.55 against predicted 1.87, 3.53, 1.89), so all three see the same underlying heterogeneity and differ only in power to resolve it. **Kv4-complex composition is not a cell-type property**: it varies more cell-to-cell within a transcriptomic type than between types, and every `*_variant_fracs_by_type.csv` table in this project reports a mean over a heterogeneous population. Within a type, the functionally decisive DPP variants are **all-or-none rather than graded** — `Vip Gpc3` carries DPP6K at a mean of 10% with 76% of cells below 5% and 8% above 50%; `L6b Serpinb11` carries DPP10a at a mean of 49% with 37% below 5% and 58% above 50% — with excess spread 30–60× the counting noise.

**A quantification artefact was found and bounded in the course of testing that claim, and it invalidates a different result.** salmon's assignment of the DPP N-terminal variants is *exact* under full-length coverage — 100% correct for all seven variants, verified at 51,000–306,000 reads and again at 385–705 reads, in both SE 100 bp and PE 50 bp — but under 3′ coverage bias it collapses into two default bins, one per gene: DPP6K and DPP6L flow into **DPP6S**, DPP10a flows into **DPP10d**. Applying the measured confusion matrix to the efficiency-calibrated qRT-PCR composition of cerebellar granule cells (Jerng & Pfaffinger 2012) predicts DPP6S at 50–57% against the 53.7% actually observed in the V404M bulk, so the artefact is real and quantitatively sufficient. Consequently **the V404M bulk DPP variant results are withdrawn**, while the single-cell results stand: SMART-seq is full-length, shows none of the collapse signature, and recovers DPP10a at 62.5% of the DPP panel in excitatory cells, consistent with DPP10a being the dominant N-type variant in cortical pyramidal neurons.

Four methodological failures found and fixed are recorded in *Technical findings*, because each produced a confident wrong result that looked like a finding. The fourth is the most instructive: the first three were caught by testing against synthetic cases, which can only catch errors already imagined, whereas this one was caught by an external measurement with a known answer.

## Datasets and methods

| Dataset | Design | Libraries | Role |
| --- | --- | --- | --- |
| V404M bulk (in-house, unpublished) | Kv4.2-V404M non-conditional vs wild-type littermate pairs, 3 months, hippocampus / cortex / thalamus / cerebellum | 48 paired-end, strand ISR (frac_fwd ≈ 0.06), mean salmon mapping 65.1% | All V404M results |
| Tasic 2016 (GSE71585) | V1, SMART-Seq v1, 100 bp SE, ~8.7 M reads/cell | 1,677 cells, 49 types | Deepest; drives the heterogeneity result |
| Allen SSv4 hippocampal formation (Yao 2021, NeMO) | SMART-Seq v4, 50 bp PE, ~0.70 M reads/cell | 8,131 cells, 26 subclasses | Heterogeneity replication at low depth |
| Allen SSv4 VISp (Tasic 2018 cells, NeMO) | SMART-Seq v4, 50 bp PE, ~2.5 M reads/cell | 1,511 cells | Intermediate depth; the √depth check |

Reference GENCODE vM36 / GRCm39 throughout; all coordinates are GRCm39. salmon 1.10.3 selective alignment, STAR 2.7.11b, featureCounts 2.1.1, samtools 1.24, pydeseq2 0.5.4. Allele counting is done on whole haplotypes straight from the FASTQs, not from alignment, because a read carrying the mutation has one extra mismatch and aligns marginally less well — an alignment-based allele ratio is biased against the mutant.

## Results, with reliability

### Allele-specific Kcnd2 expression — solid

Mutants carry the knock-in on 0.495 of Kcnd2 reads (per tissue: hippocampus 0.498, thalamus 0.504, cerebellum 0.498, cortex 0.455); wild types 0.000 in every tissue. Heterozygous, both alleles transcribed and stable in equal measure. Anchored-read recovery mutant/wild-type is 0.81, confirming both alleles are counted rather than one being missed.

This required the whole-haplotype method. The left anchor reads **0.000** in mutants because it straddles the silent AgeI site introduced beside the V404M codon; anchor-only counting would have reported the knock-in allele as transcriptionally silent. Cortex is the weakest tissue here (0.455, and anchored-read recovery 0.54 against ~0.97 in hippocampus) and its allele ratio should be treated as less certain than the other three.

### ISR, UPR and cell death — negative, but read the caveat

No panel for ATF4 output, ISR kinases, ATF6/chaperones, IRE1/XBP1, intrinsic apoptosis, p53/DNA damage, necroptosis or autophagy separates mutant from wild type in any tissue (salmon pass, pydeseq2 `~sex + genotype`, competitive and self-contained tests). The STAR junction readout agrees: XBP1 spliced-fraction medians are 0.0030 vs 0.0021 (hippocampus), 0.0059 vs 0.0049 (cortex), 0.0039 vs 0.0067 (thalamus, p = 0.022 — mutant *lower*), 0.0024 vs 0.0018 (cerebellum).

**Caveat.** The XBP1 readout rests on 220 reads pooled across all 48 libraries, with per-sample medians of 2–8 reads. That cannot support a ratio estimate. It corroborates the salmon-based negative; it does not independently establish it.

### Kcnd2 exon:intron — negative, and the sign is against the hypothesis

Intronic fraction of Kcnd2 reads, mutant vs wild type: hippocampus 0.0292 vs 0.0390 (p = 0.0221), cortex 0.0717 vs 0.0661, thalamus 0.0695 vs 0.0818, cerebellum 0.1164 vs 0.1034. The hippocampal value is the smallest p among ~36 gene×tissue tests, giving q ≈ 0.8, and the effect is one percentage point. The direction — less intron retention in mutants — is the opposite of what an increased Kv4-null fraction would produce. No population-level support for the single-cell null phenotype being elevated in V404M.

### Gene-level Kv4 changes in cortex are cell composition, not regulation — important

Cortex appears to show coordinated Kv4 downregulation: Kcnd1 −0.27, Kcnd2 −0.39, Kcnd3 −0.38, Kcnip2 −0.37, Dpp6 −0.31. It does not. Marker-class log2FC, mutant/wild type:

| tissue | neuronal | astro (identity) | microglia | oligo | Kv4 panel |
| --- | --- | --- | --- | --- | --- |
| hippocampus | −0.137 | −0.038 | +0.149 | −0.199 | +0.014 |
| cortex | **−0.185** | **+0.131** | +0.042 | +0.035 | −0.136 |
| thalamus | −0.025 | −0.025 | −0.053 | +0.017 | +0.031 |
| cerebellum | −0.102 | +0.002 | +0.009 | −0.012 | −0.086 |

The Kv4 panel moves −0.136 in cortex, *less* than the neuronal markers at −0.185, so Kv4 genes are diluted rather than selectively reduced. **Any Kv4 fold-change in cortex requires a neuronal-fraction covariate before it means anything** — including the cortex Kv4.3L PSI result (mutant 0.508 vs 0.389, p = 0.035 uncorrected), which should not be pursued as it stands. Hippocampus is unaffected by this: astrocyte identity is flat there and the Kv4 panel sits at +0.014.

### Astrogliosis: a state change, not an expansion, and probably not pathogenic — moderately solid

Reactive genes rise in hippocampus (Timp1 +0.94, Serpina3n +0.89, Gfap +0.79, Ptx3 +0.76, Hspb1 +0.65, Vim +0.43, Osmr +0.43) and more weakly in cortex (Gfap +0.62, Serpina3n +0.54, Cd44 +0.25); thalamus and cerebellum show nothing. Gfap in absolute terms: hippocampus 481 vs 313 CPM, cortex 332 vs 208.

Four observations constrain what this is.

1. **Not proliferation.** Astrocyte identity genes are flat (hippocampus −0.04, cortex +0.13) while Gfap rises 0.79. A constant astrocyte population changing state. This also rules out dissection artefact: white-matter or glial contamination would raise all astrocyte genes together and oligodendrocyte markers with them, and oligo markers instead *fall* 0.20 in hippocampus.
2. **Not the complement/MHC program.** C3 −0.24, Gbp2 and Amigo2 flat. (The A1/A2 dichotomy is best avoided altogether per the 2021 Escartin consensus; the point is that the specific neurotoxic markers are absent.)
3. **Homeostatic capacity preserved.** Per astrocyte (gene log2FC minus astrocyte-identity log2FC): Kcnj10 +0.17 hippocampus / −0.08 cortex; Slc1a2 −0.01 / −0.31; Glul −0.12 / −0.10; Gja1 −0.18 / −0.11. Kir4.1 loss (K⁺ clearance) and GLT-1 loss (glutamate clearance) are the two canonical routes by which reactive astrocytes amplify hyperexcitability, and neither is engaged. Kir4.1 is the one that would compound an A-current deficit, and it is if anything up.
4. **The regional pattern tracks circuit vulnerability, not mutant dose.** Kcnd2 median CPM in wild type: cerebellum 341, hippocampus 166, thalamus 112, cortex 101. Cerebellum expresses the most Kv4.2 by a wide margin and has no reactive signature; cortex expresses the least and does. A cell-autonomous proteotoxic response would track dose. This one tracks which circuits are epileptogenic, and is consistent with the ISR/UPR negative.

Combined with ΔFosB elevated in hippocampus (a marker of *repeated intermittent* activity, accumulating through protein half-life) and acute IEGs suppressed in all four regions, the pattern reads as a history of episodic hyperactivity sampled between episodes. That synthesis is held loosely: the ΔFosB arm is underpowered and IEG suppression could reflect handling or time of sacrifice.

**Remaining uncertainty.** Bulk RNA averages over all astrocytes, so a small severely reactive subpopulation near a focus would be invisible. And Serpina3n/Timp1/Osmr elevation is real and has been linked to synaptic remodelling, so a contribution to altered circuit function cannot be excluded even with buffering intact.

### Astrocytes do not express Kv4.2; OPCs do — a tangent worth following

Prompted by the apparent Kv4.2 preference of `Astro Aqp4` in the fraction tables, absolute levels in Tasic 2016 (reads per cell):

| type | Kv4 total | Kv4.2 | Kv4.2 share |
| --- | --- | --- | --- |
| highest neurons | 1,300–1,900 | — | — |
| **OPC Pdgfra** | **1,305.9** | **846.5** | 65% |
| SMC Myl9 | 201.3 | 193.6 | 96% |
| Micro Ctss | 168.3 | 0.3 | 0% |
| Astro Aqp4 | 11.0 | 6.5 | 59% |
| Oligo Opalin | 0.7 | 0.6 | 85% |

Astrocytes sit at 1–3% of neuronal levels; their 59% Kv4.2 *share* is a high proportion of a near-zero pool, the same compositional trap that produced the L6 CT error. **A cell-autonomous astrocytic mechanism for V404M is therefore not available** — there is not enough Kv4.2 in astrocytes for one.

OPCs are a different matter: Kv4.2 at neuronal levels, in a cell type that receives genuine synaptic input and has documented A-type currents. Microglia are a useful internal control against pervasive neuronal contamination (168 Kv4 reads/cell but essentially no Kv4.2). In the V404M data, cortical OPC markers fall coordinately (Sox10 −0.36\*, Olig2 −0.33\*, Pdgfra −0.30, Cspg4 −0.30, Gpr17 −0.27; class −0.28) while **myelinating oligodendrocytes are flat** (Mbp +0.07, Mog +0.02, Mag −0.06) — a selective precursor-pool effect, larger than the general neuronal shift. But hippocampus shows nothing (+0.08) despite having the strongest astrogliosis, so OPC and astrocyte effects dissociate regionally and this is a hypothesis, not a result. n = 22 OPCs in Tasic 2016 and a contamination check still outstanding.

### Genetic artefact controls — clean

**No passenger locus.** A knock-in targeted in 129 ES cells leaves a congenic interval around the edited gene whose genes differ by strain rather than by mutation; the signature is clustering near Kcnd2 *and* consistency across all tissues. Within ±2 Mb and ±5 Mb of Kcnd2, **zero** genes change in the same direction in ≥3 tissues (genome-wide rate 0.004), and median |mean log2FC| matches the genome-wide null at every window from 1 to 40 Mb. Kcnd2 itself: mean log2FC −0.16. The single nominal blip (±5 Mb, 4 genes at padj < 0.1, p = 0.025) fails because ±1 Mb contains zero significant genes, where an interval would be loudest.

**Sex balanced** (3F/3M mutant vs 3F/4M wild type in three tissues; cerebellum 3F/3M vs 2F/1M), with only a 1.47× chrX excess of large changes.

**The large-fold-change tail is low-count artefact.** Of 27 cross-tissue-consistent changes above |log2FC| 1, ~19 are `Gm*` loci, processed pseudogenes or bare Ensembl IDs, and most of the rest are liver/detox genes with near-zero brain expression. Gene lists from this dataset must be expression-filtered first. Two corrections to earlier verbal reports: **Gm14703 is +0.37, not −8.5** (that figure came from the salmon pass and does not survive in STAR), and **Onecut1 is hippocampus-only** (−2.50, NaN in cerebellum and cortex, +0.18 thalamus) and sits on chr9.

### Within-type heterogeneity of Kv4-complex composition — solid, and consequential

Within-subclass spread divided by the per-cell measurement noise floor, where 1.0 is pure sampling noise:

| dataset | depth/cell | Kv4 | KChIP | DPP |
| --- | --- | --- | --- | --- |
| Tasic 2016 | 8.7 M | 13.9–24.6 | 14.8–15.5 | 17.2–19.3 |
| Allen VISp | 2.5 M | 8.5–11.1 | 7.3–7.8 | 8.5–9.4 |
| Allen HPF | 0.70 M | 4.9–7.4 | 4.3–6.0 | 5.2–7.0 |

Every subclass in every dataset clears the 1.5 threshold; the most trustworthy entries are those with real cell counts (Allen HPF `DG` KChIP 5.03 at n = 1,393; `CA1-ProS` KChIP 4.13 at n = 789). The one sub-1 value is `L5/6 NP CTX` Kv4 at 0.82, n = 37.

The ratios scale as **√depth** across datasets, as they must if the biology is the same and only resolution differs: T16/VISp observed 2.02 against predicted 1.87, T16/HPF 3.13 against 3.53, VISp/HPF 1.55 against 1.89. Three datasets, two chemistries, different labs, same underlying heterogeneity.

This reframes the L6 CT problem rather than overturning it. The 20–35 percentage-point spread between `L6a Mgp`, `L6a Sla` and `L6a Car12` (Kv4.3S 3.7% / 23.6% / 38.7%; Kv4.2 81.9% / 70.6% / 60.7%) is real and our hand-made pooling into "L6 CT" was wrong. But each type individually already carries larger internal variance — `L6a Sla` alone reads 15.6 (Kv4), 10.6 (KChIP), 12.1 (DPP). Between-type structure exists; within-type scatter is bigger.

**Consequence for the project.** Kv4-complex composition is not a cell-type property. Every by-type variant-fraction table reports a mean over a population whose members differ substantially, and that should be stated wherever those tables appear. It also recasts the Kv4-null DG cells: less a distinct subpopulation than the tail of a continuum running through every subclass.

### Within a type, the decisive DPP variants are all-or-none, not graded — solid

The variance ratio above does not say what the distribution *looks like*, and the two possibilities predict entirely different electrophysiology. A population mean of 20% DPP6K is equally consistent with every cell at 20% (narrow V½-inactivation, no surprise in recordings) and with 20% of cells near 100% and the rest at zero (bimodal V½ spanning most of the DPP6 shift). Synthetic data at mean 0.20 gives excess SD 0.013 for the first case and 0.393 for the second, so the two are cleanly separable once per-cell measurement uncertainty is removed.

The data are the second case. Representative rows, where excess SD is the across-cell spread after removing counting noise:

| dataset | type | variant | mean | excess SD | noise SD | <5% | >50% | ΔBIC | shape |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| T16 | Vip Gpc3 | DPP6K | 0.103 | 0.241 | 0.001 | 0.76 | 0.08 | 245 | bimodal |
| T16 | L6b Serpinb11 | DPP10a | 0.491 | 0.419 | 0.007 | 0.37 | 0.58 | 20 | bimodal |
| T16 | L4 Ctxn3 | DPP10a | 0.397 | 0.380 | 0.006 | 0.35 | 0.42 | 198 | bimodal |
| T16 | Vip Chat | DPP6S | 0.318 | 0.339 | 0.008 | 0.30 | 0.33 | 14 | bimodal |
| T16 | Ndnf Car4 | DPP6K | 0.165 | 0.297 | 0.001 | 0.64 | 0.14 | 16 | bimodal |

Excess SD runs 30–60× the measurement noise in the deep data. Dropout cannot produce this: a cell in the >50% bin carries ≥100 reads on that one variant given the ≥200-read gate, so the split is between *absent* and *dominant*, not between absent and barely-detected. The Allen sets, with 3–4× higher noise, call the same types "graded" while still showing large excess SD — i.e. the bimodality is resolved wherever there is power to resolve it.

Predicted 1 SD spreads in V½ inactivation, under a deliberately crude linear subunit-share mapping, are 6.0–8.5 mV, implying a ±1 SD range of 12–17 mV within a single population.

### Variant assignment is exact at full-length coverage and collapses into two sinks under 3′ bias — solid, and it invalidates the bulk DPP results

Nothing in this project had ever checked salmon's isoform assignment against a known ratio. The 2012 standardizing constructs — tandem plasmids carrying GAPDH exons 5–8 followed by exons 1–3 of two DPLP variants at a guaranteed 1:1:1 molar ratio — are precisely that control for qRT-PCR, and the logic transfers: simulate equal molar amounts of each variant, quantify against this project's own index, and compare.

**At full-length coverage, assignment is exact**: 100% correct for all seven DPP variants, at 51,000–306,000 reads per variant and again at 385–705 reads, in SE 100 bp (the Tasic 2016 layout) and PE 50 bp (the Allen layout). There is no intrinsic sequence ambiguity and no low-depth instability.

**Under 3′ coverage bias it collapses into two default bins, one per gene** (PE, tau = 250 nt):

| true origin | correctly assigned | where the rest goes |
| --- | --- | --- |
| DPP6K | 0–34% | **66–100% → DPP6S** |
| DPP6L | 37% | **63% → DPP6S** |
| DPP10a | 23% | **77% → DPP10d** |
| DPP6S, DPP6a, DPP10c, DPP10d | 100% | — |

DPP6S absorbs 129–166% of other variants' reads summed; DPP10d absorbs 77% of DPP10a's. The magnitude is sufficient to explain the real data: pushing the 2012 cerebellar granule cell composition (DPP6K 41%, DPP6a 33%, DPP6S 18%, DPP6L 8%) through this matrix predicts DPP6S at 50–57% against 53.7% observed in the V404M bulk cerebellum, where an unbiased pipeline would have returned 18%.

**Therefore the V404M bulk DPP variant results are withdrawn.** The bulk libraries are 3′-biased polyA, and they report DPP10a at 12.4% of the DPP panel in cortex where the full-length single-cell data gives 62.5% — a ~5-fold under-call of the variant that carries N-type inactivation in cortical pyramidal neurons. DPP6K is destroyed outright. The one bulk DPP measurement that survives is **DPP6a**, which never collapses at any tau; its near-absence in cortex (0.44% of the panel) is therefore real, and consistent with DPP10a rather than DPP6a carrying the N-type function there.

The single-cell results are unaffected and show none of the collapse signature — DPP10a reads 62.5% in excitatory cells with DPP10d at 15.4%, the opposite of the collapse direction:

| | DPP10a | DPP10d | DPP6S | DPP6a | DPP6K |
| --- | --- | --- | --- | --- | --- |
| Tasic 2016 excitatory | 62.5% | 15.4 | 8.6 | 8.3 | 0.0 |
| Tasic 2016 inhibitory | 2.8 | 35.8 | 12.7 | 31.0 | 9.7 |

One route remained by which the artefact could have reached the single-cell data: RNA integrity varies cell to cell, so an intact cell reports DPP6K correctly while a degraded one reports the same molecule as DPP6S, splitting cells into "expresses it" and "doesn't" by quality alone — which is the observed bimodality. Because Dpp6 and Dpp10 have *independent* sinks, a shared coverage profile would move both together. It does not: across 29 cell types in three datasets, the DPP6S–DPP10d correlation has a median rho of −0.005 and sits a median 0.254 *below* its compositional null, with 1 of 29 types exceeding the null (chance alone gives ~5%). A synthetic shared-quality driver produces +0.63 beyond null, so this is not a power limitation. **The bimodality is not a coverage artefact.**

### What a pipette would encounter, and what it would take to test it — solid

Recast as recordable categories, with the A-current-bearing population (any Kv4 detected) as the denominator,
since that is the population a patch electrode samples from. For `L6a Sla` (n = 62):

| | % |
| --- | --- |
| cells with any Kv4 subunit | 93.5 |
| of those: full ternary complex (KChIP + DPP) | 96.6 |
| of those: DPP only, no KChIP | 1.7 |
| of those: KChIP only, no DPP | 1.7 |
| **of those: no Kv4.2 (Kv4.1/Kv4.3 only)** | **20.7** |
| of the Kv4.2-negative cells, those carrying the full auxiliary set | 91.7 |

**The auxiliary complement does not vary in cortical pyramidal cells.** Across every excitatory type the
A-current-bearing cells carry KChIP + DPP essentially without exception — L2/3 Ptgs2 98.9%, L4 Arf5 100%,
L4 Scnn1a 100%, L5a Batf3 100%, L5b Tph2 100%, L6a Mgp 100%, L6a Car12 100%, L6b Serpinb11 100%. What varies
is the pore subunit, not its partners. An earlier version of this analysis reported ~73% "full complex" for
L6a Sla; that was an artefact of a category scheme in which categories 3–6 all required Kcnd2 and so never
assessed the auxiliary complement of a cell running on Kv4.1/Kv4.3. On synthetic data where 97.7% of
A-current cells carried the full complement, that scheme returned 73.6%.

**Partial complexes do exist, but not in cortex.** `77_Sst HPF` has 73.3% of its A-current cells on DPP
without KChIP, `78_Sst HPF` 37.0%, `8_Lamp5 Lhx6` 36.0% — a real hippocampal interneuron configuration, with
the faster kinetics and hyperpolarised V½ that implies.

**Within L6a, Sla carries the largest Kv4.2-negative fraction** (20.7%, against Mgp 11.1% and Car12 5.0%).
Note the Allen data puts `L6 CT CTX` at 3.7–6.2% and `L6 IT CTX` at 16.7%, so if `L6a Sla` is IT rather than
CT — and the Oct 30 mapping of L6a types onto CT/IT was our inference, not a published crosswalk — the two
datasets reconcile.

**Implication for V404M.** About 19% of L6a Sla cells one would patch carry no Kv4.2 and should be
genotype-insensitive: identical between V404M and wild-type littermates. That predicts a **mixture** in the
mutant A-current distribution against a unimodal wild-type one — a within-animal control needing no
sequencing.

**How many recordings.** Simulated, using wild-type littermates as a fixed reference component and testing by
likelihood ratio whether the mutant sample needs it (far more powerful than asking one sample whether it
looks bimodal; critical values calibrated by simulation because the mixing proportion is on a boundary under
H0). With d = shift ÷ within-type SD of the parameter:

| | power at n = 80 | power at n = 200 | estimate of p at n = 200 |
| --- | --- | --- | --- |
| d = 1.0 | 0.06 | 0.13 | — |
| d = 1.5 | 0.21 | 0.44 | 0.19 [0.00–0.32] |
| d = 2.0 | 0.60 | 0.76 | 0.19 [0.08–0.28] |
| d = 3.0 | 0.97 | 1.00 | 0.19 [0.14–0.25] |

**The hypothesis test is the wrong goal; the estimate is the right one.** Formal detection of the mixture
needs d ≥ 3 (80% power at n = 40) or n beyond 200 at d = 2. But the *estimate* of the unshifted fraction is
unbiased from n ≈ 80 — median 0.19 against a true 0.194 at every effect size — and its precision is what
decides feasibility: at d = 3, n = 80 gives 0.19 [0.10–0.29]; at d = 2, n = 200 gives [0.08–0.28]; at
d = 1.5 the interval includes zero even at n = 300, so the experiment is not worth attempting.

**Therefore the first step is to measure d, not to plan the experiment.** The lab's existing V404M and
wild-type L6 recordings give it directly: the genotype difference in mean V½ inactivation divided by the
within-type SD. That single number says whether this is an 80-cell experiment, a 200-cell experiment, or one
not to run.

**What fraction truly lacks an A-current — resolved, and it is essentially none.** Among the deeper half of
L6a Sla libraries (n = 31), **not one cell has zero Kv4 reads** (95% upper bound 9.2%), against 6.5% at the
≥5-read threshold across all cells. The 6.5% was detection. The metric validates on the controls, which stay
zero at depth where they should: Astro Aqp4 63.6%, Oligo Opalin 53.3%, Micro Ctss 54.5%.

**But the Kv4.2-negative population survives the same test, and that dissociation is the key result.** In the
identical cells, 19.4% (6 of 31) carry zero Kcnd2 reads while every one of them has Kv4 reads from Kcnd1 or
Kcnd3 — against 20.7% at the ≥5-read threshold, i.e. essentially unchanged. One call collapsed with depth
and the other did not, which is strong evidence the second is biological: these are cells running an
A-current on a different pore subunit, not cells sequenced too shallowly. `L5a Batf3` shows the same
dissociation (0.0% Kv4-null, 21.2% Kv4.2-null) and `L6a Mgp` likewise (4.8% / 14.3%).

It is type-specific rather than general: `L5b Tph2`, `L5b Cdh13`, `L5a Tcerg1l` and `L6a Car12` sit at 0%,
though with n of 10–18 their upper bounds are 15–26% and a modest population is not excluded. Among
interneurons the same test confirms **`Pvalb Tacr3` at 5.6%** (UB 16.5%) and `Sst Chodl` at 4.8% as the
cleanest targets, against `Pvalb Rspo2` 81.2% and `Vip Gpc3` 80.0%.

**So the corrected figures for an L6a Sla recording are: ~0% with no A-current, ~19% with an A-current that
V404M cannot touch, ~81% affected.** That is the p = 0.194 the power analysis assumed, now established
rather than asserted. The earlier statement that 6.5% should have no A-current is superseded.

The residual caveat is unchanged and unquantified: bursting means a cell with no detectable Kcnd2 mRNA may
still hold Kv4.2 protein made earlier, so 19% remains an upper bound on the genuinely genotype-insensitive
fraction, by an amount set by the protein half-life relative to the message.

*(Superseded text follows for the record:)* 6.5% was an upper bound. Detection
threshold, transcriptional bursting and an arbitrary cutoff all push the same way, and mature Kcnd2 mRNA is a
long-lived pool fed by bursts (companion document, Oct 30), so a cell silent at capture may still hold
channel protein. The defensible estimate is the fraction with exactly zero Kv4 reads among the deepest
libraries — a cell sequenced to millions of reads with not one Kv4 read — computed by
`scripts/recording_design.py`.

### Effect dilution: why a null result in some interneurons would be uninterpretable — important for design

Kv4.2 is far less dominant in interneurons, so a larger share of any interneuron sample is genotype-
insensitive. That dilutes the apparent effect twice over: the mean shift falls to (1−p)·Δ while the SD rises
to √(σ² + p(1−p)Δ²), since the mixture adds between-component variance. Taking a true per-cell effect of
d = 2:

| recorded type | % no Kv4.2 | apparent d | n per group, 80% power |
| --- | --- | --- | --- |
| **Pvalb Tacr3** | **1.6** | **1.91** | **5** |
| Sst Chodl | 12.2 | 1.47 | 8 |
| L5a Batf3 | 19.3 | 1.27 | 10 |
| L6a Sla | 20.7 | 1.23 | 11 |
| Sst Myh8 | 33.3 | 0.97 | 17 |
| Vip Parm1 | 44.1 | 0.79 | 25 |
| Vip Gpc3 | 72.7 | 0.41 | 95 |
| Lamp5 Lhx6 (HPF) | 82.9 | 0.27 | 211 |
| Vip HPF | 92.0 | 0.14 | 793 |

Three consequences.

**A null result in Vip, Sncg or hippocampal Lamp5 types carries almost no information** at ordinary sample
sizes: a large per-cell effect and no effect at all look the same with 10–20 cells per group. Most
hippocampal interneuron classes are in this regime.

**Interneurons are not inherently poor targets — the problem is type-specific.** `Pvalb Tacr3` at 1.6%
Kv4.2-negative would need five cells per group, fewer than any pyramidal type, and `Sst Chodl` eight. If the
question is whether V404M affects inhibition, those are the recording targets.

**The bias has a direction, so existing data is readable.** Dilution can only shrink an effect, never create
one, so any significant interneuron effect already in hand is an **underestimate of the per-cell effect by
1/(1−p)** — 3.7× for Vip Gpc3, 5.8× for hippocampal Lamp5 Lhx6.

The remedy that removes dilution entirely is to genotype the recorded cell: single-cell RT-PCR for Kcnd2 on
the harvested cytoplasm, analysing only the cells the mutation can act on. The 2012 nested primer design
does exactly this, and needs only presence/absence of one transcript.

**Cre lines.** There is no Tacr3-Cre, and cluster-level driver lines generally do not exist below the
subclass. The practical options are the standard subclass lines — Pvalb-IRES-Cre, Sst-IRES-Cre,
Vip-IRES-Cre, Ndnf-IRES2-dgCre — followed by post-hoc identification, or an intersectional Cre × Flp
strategy. Note that PV-IRES-Cre has documented incomplete labelling of Pvalb cells in some cortical regions,
so labelling efficiency should be checked for the region of interest rather than assumed. The empirical
question of which lines actually yield which cluster is answerable from the Tasic 2016 metadata, since every
cell was FACS-sorted from a named driver line crossed to tdTomato. That mapping (1,530 of 1,677 cells):

| target | best line | yield | cells to patch per 5 of target |
| --- | --- | --- | --- |
| Sst Chodl (4.8% Kv4.2-neg) | **Nos1** | 43.1% (31/72) | ~12 |
| Pvalb Tacr3 (5.6%) | PvalbF-Gad2 | 32.2% (28/87) | ~16 |
| Pvalb Tacr3 | **Pvalb** | 27.0% (24/89) | ~19 |
| L6a Sla (19.4%) | **Ntsr1** | 61.0% (47/77) | ~8 |
| L2/3 Ptgs2 | Cux2 | 47.5% (58/122) | ~11 |
| L5a Batf3 | Rbp4 | 23.1% (40/173) | ~22 |

Plain Pvalb-IRES-Cre at 27% makes the intersectional `PvalbF-Gad2` (32%) not worth the extra cross, and
`Nos1` at 43% is simply hitting its native target, since `Sst Chodl` is the long-range-projecting Sst/NOS1
population. **The dilution-aware interneuron experiment is therefore cheap**: at d = 2, `Pvalb Tacr3` needs
five cells per group, or about 19 Pvalb⁺ cells patched per genotype.

Yields are rank orderings, not hit rates — Tasic 2016 stratified sampling by line and dissected layer — and
they come from VISp, so composition may differ in other cortical areas.

**One consequence worth recording.** `Ntsr1` is the canonical L6 corticothalamic driver and supplies 78.3%
of `L6a Sla` cells, so **L6a Sla is CT**. The earlier suggestion that the gap between `L6a Sla` (19.4%
Kv4.2-negative) and Allen `L6 CT CTX` (3.7–6.2%) might reconcile if L6a Sla were IT therefore fails. The
discrepancy is real and unexplained; depth does not account for it either, since the Allen data is shallower
and should over-call rather than under-call Kv4.2-negative. The likeliest remaining explanation is that the
two clusters are not the same cells despite both being broadly CT — the taxonomy mismatch recorded above.

### Technical findings worth recording

**Four methodological failures, each of which produced a confident wrong answer.**

*Centred log-ratio with a pseudocount reads out sequencing depth on sparse data.* A structural zero contributes log(pseudo/total), so two cells with identical composition at 100 vs 5,000 reads differ by log(50) ≈ 3.9 in every zero entry. Synthetic cells with identical composition, true zeros and depths spread over 100–5,000 returned a spread/noise ratio of **12.23 against a true value of 1.0**; the same cells gave 0.85 under Hellinger coordinates (√proportions, where √0 = 0 at any depth). Sparsity and depth variation are both required to trigger it, which is exactly single-cell isoform data. This invalidated the first grouping run, which returned "split" for 24 of 28 groupings — including Allen subclasses against their own clusters — on R² of 0.01–0.22.

*Split-half resampling cannot estimate noise for concentrated compositions.* It conditions on the observed support, so a cell with all reads in one variant yields **identically zero** noise at any depth (normalising [n,0,0,0] always gives [1,0,0,0]). Real panels are full of such cells; each contributed a zero to the median and the ratio inflated without bound (observed 8–211, where Hellinger distance to a centroid is bounded at 1.0, so those values were arithmetically impossible). Replaced with a Dirichlet(counts + 0.5) posterior, which represents uncertainty about *unobserved* variants: for [200,0,0,0] the truth could be [0.995,0.005,0,0]. Validated at 1.02 on diffuse and on concentrated Kv4-like nulls at variable depth, conservative (0.34) on very sparse compositions, and sensitive to real heterogeneity (1.55/3.34 mild, 3.25/8.85 strong at depth 200/2000).

*Mutual exclusivity cannot distinguish real exclusive expression from EM collapse.* A test built to separate them could not, because exclusivity is the expected signature of both: a cell genuinely expressing only DPP6K is indistinguishable from one whose reads the EM dumped there. Its validation hid this by constructing the "biological" case as cells in one of two *mixture* states, so the answer was built into the premise, and it also flagged its own Kv4.3S/L positive control. Abandoned.

*Resampling measures precision, not accuracy, and the error here was accuracy.* Bootstraps and split-half both resample reads **under the quantifier's own model**, so a systematic model failure is invisible to them. The planned bootstrap run — two machines overnight, ~90 GB of re-downloaded FASTQ each — would have returned a reassuring answer while the DPP6S bias went undetected, exactly as the 4-fold Kv4.3L underestimate from 3′UTR model asymmetry would not appear in a bootstrap. It was cancelled in favour of the known-ratio calibration, which cost an hour.

*General lesson, and the one worth carrying forward.* Compositional methods needing a pseudocount are hazardous on sparse, depth-variable data, and any new metric should carry a known-answer null — the noise-floor metric's failures surfaced only because pure noise *must* give 1.0, making 131 visibly impossible, whereas the original R² had no such anchor and its failure looked like a result. But note the asymmetry between the first three failures and the fourth. The first three were caught by testing against synthetic cases, which can only catch errors already imagined. The DPP assignment bias was caught by an **external measurement with a known answer** — the lab's own 2012 efficiency-calibrated qRT-PCR — and no amount of internal self-consistency checking would have found it, because every internal check shares the quantifier's assumptions. Isoform-level claims in this project should be anchored to an outside standard before they are believed.

**Infrastructure.** `hip_common.py` pointed `HIP_PKL`/`CTX_PKL` at paths inside an ephemeral analysis container; those pickles no longer exist anywhere. Both loaders now read the durable annotated CSVs from the collapse step (`HIP_DATA`/`CTX_DATA`). `region_compare3.py` still carries a dead scratchpad path and needs the same treatment before it can be re-run.

## Is the within-type heterogeneity biological or methodological?

Largely settled, in favour of biological, by four tests that between them close the routes a technical
explanation could take. The honest residue is stated at the end.

**Sampling noise is excluded.** The per-cell Dirichlet posterior puts the within-subclass spread at 4–25×
the counting noise, and the ratio scales as √depth across three datasets and two chemistries, so all three
see the same heterogeneity and differ only in resolution.

**Isoform-assignment error is excluded for this panel, by direct calibration.** This was the main worry, and
it was well-founded — the artefact is real and large, but it requires 3′ coverage bias. Under full-length
coverage assignment is exact (100% for all seven DPP variants) at both the high depths and the single-cell
depths where the heterogeneity is measured, and in both the SE 100 bp and PE 50 bp layouts the datasets
actually use. The single-cell data show none of the collapse signature.

**Cell-to-cell variation in RNA integrity is excluded.** This was the one way the assignment artefact could
still have reached the single-cell data, and it is the mechanism that would most naturally manufacture
all-or-none variant usage. Because Dpp6 and Dpp10 have independent sinks, a shared per-cell coverage profile
would move both together; across 29 cell types the sink–sink correlation sits a median 0.254 *below* its
compositional null, where a synthetic shared driver gives +0.63 above it.

**Dropout is excluded by the shape of the distribution.** Cells in the >50% bin carry ≥100 reads on that one
variant, so the split is between absent and dominant, not between absent and barely-detected.

What remains unexcluded, in rough order of concern:

1. **PCR amplification noise.** SMART-seq amplifies before sequencing, and stochastic early-cycle
   amplification is multiplicative and depth-independent, so it sits in the numerator and the noise floor
   cannot see it. Both v1 and v4 amplify, so cross-chemistry agreement does not rule it out. A negative-
   control gene family would address this: run the identical analysis on isoform families of comparable
   expression with no expected biological heterogeneity. If those return ratios of 10–20 too, the variance
   is generic and technical; if the Kv4-complex ratios are specifically high, that is strong evidence for
   biology. Uses existing data, a few hours, and is the one remaining cheap test.
2. **Reverse-transcription efficiency** differing between variants and varying with cell state.
3. **Dissociation-induced transcriptional change**, which both protocols impose and which can shift splicing.

And biological itself splits two ways with different consequences: **stable cell-to-cell differences**
(developmental, epigenetic, allelic) versus **dynamic state** (activity, circadian, cell cycle). The first
makes a cell's A-current phenotype a durable property; the second makes it fluctuate within a cell over
hours. Nothing here distinguishes them, and the distinction matters for whether composition could be
reprogrammed.

The defensible statement is: *within-type variance exceeds sampling noise by 4–25×, is not isoform-assignment
error, is not coverage-driven, and takes an all-or-none rather than graded form in the functionally decisive
DPP variants. Amplification noise remains formally unexcluded and is testable with existing data.*

## Predicted consequences for recorded A-currents in pyramidal neurons

The Kv4 complex is a three-component machine, and the auxiliary subunits — not the pore — set most of the biophysics. So composition variance maps onto specific, measurable current parameters rather than onto current amplitude alone. Since the composition data now survive calibration, these are predictions to test rather than conditionals.

The quantitative anchor: the measured all-or-none distributions give predicted 1 SD spreads of 6.0–8.5 mV in V½ inactivation under a linear subunit-share mapping, i.e. a ±1 SD range of 12–17 mV within one pyramidal population — large enough to dominate most other sources of variance in a recording.

**What each component contributes.** KChIPs increase surface expression, slow inactivation and speed recovery, with **KChIP4a the outlier**: its N-terminal KIS domain suppresses inactivation outright, converting fast-inactivating A-type current toward a slowly- or non-inactivating form. DPP6 and DPP10 shift activation and steady-state inactivation to markedly hyperpolarised potentials (on the order of 20–25 mV for DPP6) and accelerate both inactivation and recovery; their N-terminal variants (DPP6-S vs DPP6-L vs DPP6-K; DPP10a) differ in how strongly they do so. Kv4.2 and Kv4.3 differ more modestly from each other, and the Kv4.3L cassette adds a C-terminal insert.

**Therefore the predicted phenotype, in order of expected variance:**

1. **V½ of steady-state inactivation should be the most variable parameter**, driven by cell-to-cell differences in DPP6/DPP10 content and N-terminal variant usage. This matters disproportionately because it is measured near rest: a 10–15 mV spread in V½inact translates into a large spread in how much A-current is actually *available* at −70 mV, even between cells with identical peak conductance.
2. **Recovery from inactivation (τrec) should vary substantially**, being both KChIP- and DPPL-dependent. The functional readout is behaviour during trains — spike-frequency adaptation, dendritic summation of back-propagating action potentials, and the degree to which A-current limits repetitive firing.
3. **Inactivation kinetics and completeness should vary**, with a subset of cells showing slowed or incomplete inactivation wherever KChIP4a is present. If KChIP4a usage is bimodal rather than graded, this predicts a **discrete subpopulation** with qualitatively different current, not a smooth spread.
4. **Peak A-current density should vary**, but is the *least* diagnostic parameter, since it also depends on total channel number, trafficking and recording conditions.

**Testable signatures in existing and new recordings:**

- Within a single pyramidal type, the variance in V½inact and τrec should **exceed what recording noise and series-resistance error can account for** — an analysis that can be done retrospectively on recordings the lab already holds, with no new experiments.
- Parameters should **co-vary within cells** rather than vary independently, if composition is coherent (a cell rich in DPP6 and KChIP1 should sit at a predictable combination of V½inact and τrec). Independent scatter across parameters would argue for measurement noise; correlated structure argues for composition.
- Cluster structure in the ephys parameter space should **not respect transcriptomic type boundaries**, mirroring the RNA result.
- **Patch-seq is not the right instrument**, and the calibration says why. If the distribution is all-or-none, as it is, the subpopulations are separable by patch-clamp alone and no sequencing is needed; patch-seq would only be necessary in the graded case, where per-cell composition must be correlated against per-cell kinetics — and that is precisely where its depth (0.5–2 M reads, against 8.7 M for Tasic 2016, with RNA degraded by dialysis and partial aspiration) is least adequate. Patch-seq is most useful in the regime where it is least feasible. For presence/absence of a handful of variants, which is all a bimodality claim requires, targeted single-cell RT-PCR with the 2012 nested primers is more sensitive and involves no assignment inference at all.

**Two caveats that could still break the chain, both downstream of the RNA.** mRNA composition is not protein composition: translation efficiency, trafficking, assembly stoichiometry and protein half-life all intervene, and the lab's own evidence for direct FMRP–Kv4 protein interaction (Zhan 2024) is a reminder that regulation sits downstream of transcript abundance. And the complex assembles somatically (companion document, Oct 30), so a cell's surface channel population integrates its transcript composition over the mRNA and protein lifetime — which would *smooth* transcript-level variance and predict less current variance than the RNA suggests. If recorded A-current parameters turn out to be *tighter* than composition variance predicts, that difference is itself informative: it would imply a post-transcriptional buffering mechanism setting A-current to a target.

**Relevance to V404M.** Since auxiliary-subunit composition does vary this much between neighbouring pyramidal cells, the functional impact of V404M should vary with it — the mutation acts on a channel whose properties are set by partners the cell chooses. That predicts a *distribution* of severity across cells within a region rather than a uniform deficit, and offers a mechanism for why a single V404M copy is worse than one or two nulls: a null removes a channel, whereas V404M inserts a defective subunit into every complex the cell assembles, whatever its composition.

## Concerns and open questions

1. The cortex composition confound invalidates all cortical gene-level Kv4 calls until a neuronal-fraction covariate is added. Hippocampus is unaffected.
2. The XBP1 and ΔFosB junction readouts are too shallow to support ratio estimates and should not be quoted as quantitative.
3. PCR amplification noise is the one technical explanation for the within-type heterogeneity still formally unexcluded; a negative-control gene family would close it with existing data.
3b. **All DPP variant numbers from the V404M bulk are withdrawn** (3′-bias collapse into DPP6S and DPP10d). DPP6a is the exception and survives. Any earlier document or table quoting bulk DPP6/DPP10 composition needs the same correction.
4. The OPC Kv4.2 observation rests on 22 cells and lacks a contamination check.
5. ΔFosB cannot be assayed properly by bulk junction counting; the IHC remains the better evidence.
6. Bulk RNA cannot exclude a small, severely reactive astrocyte subpopulation.
7. `region_compare3.py` and the three-way comparison it produced are not currently reproducible (dead paths, and the CLR metric issue applies to its distance calculations).
8. The DPP calibration covers the DPP panel only. The Kv4 and KChIP panels have not been calibrated, and the Kv4.3L 3′UTR underestimate found earlier shows that panel is not immune. The same script covers them with one flag change.
9. Whether the per-cell heterogeneity is a stable property or a dynamic state is untested, and it decides whether composition could be reprogrammed.
10. The ~19% Kv4.2-negative figure assumes a cell with no detectable Kcnd2 mRNA has no Kv4.2 protein. Bursting makes that an upper bound, by an amount set by the protein half-life relative to the message, which nothing here measures.
11. The mixture experiment is infeasible below d = 2 regardless of n, so measuring d comes before designing it.

## Recommended additional analyses

| Priority | Analysis | Why |
| --- | --- | --- |
| 1 | Measure d (genotype shift ÷ within-type SD) from existing V404M and wild-type L6 recordings | One number, from data already in hand, that decides whether the mixture experiment is 80 cells, 200 cells, or not worth running |
| 1a | Re-read any existing V404M interneuron recordings in light of the dilution table | A non-significant result in Vip or Lamp5 types means little; a significant one is an underestimate by 1/(1−p) |
| 1b | Retrospective variance analysis of existing pyramidal A-current recordings | Tests the central prediction at no cost. V½ inactivation and recovery kinetics should vary beyond recording noise, and co-vary within cells |
| 2 | Negative-control gene family through the noise-floor analysis | Closes the last technical explanation for the heterogeneity; existing data, a few hours |
| 3 | Calibrate the Kv4 and KChIP panels (`simulate_variant_calibration.py --panels Kv4,KChIP`) | The DPP panel turned out to have a 5-fold bias; these are unchecked, and Kv4.3L is already known to be mis-estimated |
| 3b | Neuronal-fraction covariate / deconvolution for cortex and hippocampus | Makes the V404M gene-level numbers interpretable |
| 4 | Neuronal-specific conditional V404M | Separates neuronal from OPC contributions *and* tests whether astrogliosis is secondary — one experiment, two questions |
| 5 | Developmental time course of the astrogliosis | Cause vs consequence |
| 6 | Neuronal-marker contamination check on Tasic 2016 OPC libraries | Cheap; closes or opens the OPC thread |
| 7 | Pdgfra⁺/Olig2⁺ counts in V404M cortex | Tests the OPC prediction directly |
| 8 | Astrocyte-resolved (sorted or snRNA) sequencing with per-cell Kcnj10 | Excludes a severely reactive subpopulation |
| 9 | Rebuild `channel_phenotype_groups.py` with the Hellinger metric and the noise floor as effect-size gate | Its first run is retracted; now worth doing |
| 10 | Targeted single-cell RT-PCR for DPP variants using the 2012 nested primers | For presence/absence — which is all bimodality requires — this beats whole-transcriptome sequencing on both sensitivity and freedom from assignment inference |

### Going forward

The heterogeneity question has moved from open to largely settled, so the next step is no longer a
sequencing experiment but an electrophysiological one, and it costs nothing. If within-type DPP composition
is all-or-none as the data say, then V½ inactivation within a single pyramidal population should be spread
over roughly 12–17 mV and possibly bimodal, parameters should co-vary within cells, and the clustering
should ignore transcriptomic type boundaries. All three are testable retrospectively on recordings already
in the lab's files. The asymmetry is worth exploiting: a tight unimodal V½ inactivation in a population the
RNA calls bimodal would **falsify** the functional consequence outright, far more cheaply than any
sequencing could confirm it.

Patch-seq is not the way to pursue this, and the calibration explains why. If the distribution is
all-or-none, the subpopulations are findable by patch-clamp alone and no sequencing is needed; if it were
graded, per-cell composition would have to be correlated against per-cell kinetics at depths where DPP
recovery is marginal. Patch-seq is most useful in exactly the regime where it is least feasible. Targeted
single-cell RT-PCR with the 2012 nested primers is the better instrument for presence/absence, which is all
a bimodality claim requires.

For V404M, the biology now points away from proteotoxicity and toward circuit dysfunction: both alleles expressed, no stress response, no change in Kcnd2 processing, astrogliosis confined to epileptogenic regions and not tracking mutant dose, homeostatic astrocyte function intact. The conditional line is the experiment that would most efficiently advance this, since it addresses the neuronal-vs-OPC and the cause-vs-consequence questions at once.

## Data, code and provenance

Repository `ppfaff/Kv4-isoform-survey`. V404M pipeline: `pipeline/build_v404m_runlist.py`, `find_v404m_codon.py`, `set_v404m_haplotypes.py`, `v404m_allele_extract.py`, `run_v404m_salmon.sh`, `run_v404m_star.sh`, `collapse_v404m.py`; analysis `scripts/v404m_panels.py`, `v404m_stress.py`, `v404m_kv4.py`, `v404m_passenger.py`. Recording-facing: `scripts/complex_composition_per_cell.py` (recordable categories, pore x auxiliary cross-tab, N-type, DPP6K, detection checks), `scripts/recording_design.py` (true-absence estimate; mixture power by EM). Heterogeneity and calibration: `scripts/channel_noise_floor.py` (Dirichlet floor), `variant_distribution_shape.py` (graded vs all-or-none), `sink_coupling.py` (RNA-integrity control), `pipeline/simulate_variant_calibration.py` (known-ratio calibration, mix and purity modes), `channel_phenotype_groups.py` (grouping, retracted pending rebuild), `hip_common.py`. Abandoned: `variant_exclusivity.py` (cannot separate the hypotheses it was built for; see *Technical findings*). Not run: `pipeline/run_bootstrap_uncertainty.sh` and `scripts/bootstrap_vs_biological.py`, superseded by the known-ratio calibration before execution. Dataset provenance for all public data in `lab_db/DATASETS.md`. Output tables in `C:\Users\paulp\kv4_v404m_tables` and `tables/` in the repository.

Reference measurement: Jerng HH, Pfaffinger PJ (2012) *Incorporation of DPP6a and DPP6K variants in ternary Kv4 channel complex reconstitutes properties of A-type K current in rat cerebellar granule cells.* PLoS ONE 7(6):e38205 — the standardizing constructs described there (GAPDH exons 5–8 in tandem with exons 1–3 of two DPLP variants, 1:1:1 molar) are the only known-ratio standard this project has, and are what identified the DPP assignment bias.

Companion documents: *Kv4 Complex mRNAs Are Somatic* (Oct 30) and the isoform survey conclusions (Oct 29).
