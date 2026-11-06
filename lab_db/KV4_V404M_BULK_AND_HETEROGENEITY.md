# Kv4.2-V404M Bulk RNA-seq, and Within-Type Heterogeneity of Kv4-Complex Composition

Nov 6, 2026 · @Paul

## Summary

Two separate pieces of work, reported together because the second changes how the first should be read.

**V404M bulk RNA-seq** (48 libraries, littermate pairs, 4 tissues, 3-month-old mice) answers the four questions put to it, three of them negatively. Mutant animals are heterozygous with **both alleles expressed at parity** (knock-in fraction 0.495; the two haplotype windows carry 5,510 and 5,485 reads), so the missense change does not destabilise its own transcript. There is **no integrated stress response, no unfolded-protein response and no apoptotic program** in any tissue: the XBP1/IRE1 splice shows no mutant elevation anywhere, and the one nominal hit (thalamus, p = 0.022) runs the wrong way. The **Kcnd2 exonic:intronic ratio is unchanged**, and the direction in hippocampus is *less* intron retention in mutants (0.0292 vs 0.0390, p = 0.0221 nominal, q ≈ 0.8 across ~36 gene×tissue tests) — the opposite of an increased unprocessed/null fraction. The **ΔFosB junction** is elevated 2.3-fold in mutant hippocampus (0.0667 vs 0.0294, p = 0.100), the largest of the four tissues, with cortex flat and cerebellum at zero; the tissue-specificity matches the lab's ΔFosB immunostaining, but at 109 reads pooled across 48 libraries this is underpowered, and ΔFosB accumulates through protein half-life rather than transcript abundance, so bulk junction counts are the wrong assay for it either way.

What the mutants *do* show is **a reactive astrocyte state in hippocampus and cortex and nowhere else** (Gfap +0.79 and +0.62; Serpina3n +0.89 and +0.54; Timp1, Ptx3, Vim, Osmr in hippocampus). Three findings constrain its interpretation. It is **not proliferation**: astrocyte identity genes (Aqp4, Slc1a3, Gja1, Aldoc) are flat at −0.04 in hippocampus. It is **not the complement/neurotoxic program**: C3 is *down*. And **per-astrocyte homeostatic capacity is preserved** — Kcnj10 (Kir4.1) +0.17 in hippocampus, Slc1a2 (GLT-1) −0.01, so the two canonical routes by which reactive astrocytes amplify hyperexcitability are not engaged. Decisively, the regional pattern does **not** track mutant protein dose: cerebellum expresses the most Kcnd2 of any tissue (341 CPM vs hippocampus 166, thalamus 112, cortex 101) and shows no reactive signature at all. It tracks which circuits are epileptogenic. On this evidence the astrogliosis is a **readout of circuit dysfunction, not a driver of it**.

**Within-type heterogeneity.** Asked whether the large per-cell variance in Kv4/KChIP/DPP composition reflects lumping distinct neuron types or intrinsic variation within coherent ones, the answer is unambiguous: **intrinsic**. Every subclass in all three single-cell datasets shows within-subclass spread 4–25× the measurement noise floor. The ratios scale as √sequencing-depth across the three datasets (observed 2.02, 3.13, 1.55 against predicted 1.87, 3.53, 1.89), meaning all three see the same underlying heterogeneity and differ only in their power to resolve it. **Kv4-complex composition is therefore not a cell-type property**: it varies more cell-to-cell within a transcriptomic type than between types. Every `*_variant_fracs_by_type.csv` table in this project reports a mean over a highly heterogeneous population. Whether that heterogeneity is biological or a residue of methodology is not yet settled — sampling noise is excluded, isoform-assignment uncertainty is not — and the tests that would settle it are specified below.

Two metric failures found and fixed along the way are recorded in *Technical findings*, because both produced confident, wrong results that looked like findings.

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

### Technical findings worth recording

**Two metric failures, both of which produced confident wrong answers.**

*Centred log-ratio with a pseudocount reads out sequencing depth on sparse data.* A structural zero contributes log(pseudo/total), so two cells with identical composition at 100 vs 5,000 reads differ by log(50) ≈ 3.9 in every zero entry. Synthetic cells with identical composition, true zeros and depths spread over 100–5,000 returned a spread/noise ratio of **12.23 against a true value of 1.0**; the same cells gave 0.85 under Hellinger coordinates (√proportions, where √0 = 0 at any depth). Sparsity and depth variation are both required to trigger it, which is exactly single-cell isoform data. This invalidated the first grouping run, which returned "split" for 24 of 28 groupings — including Allen subclasses against their own clusters — on R² of 0.01–0.22.

*Split-half resampling cannot estimate noise for concentrated compositions.* It conditions on the observed support, so a cell with all reads in one variant yields **identically zero** noise at any depth (normalising [n,0,0,0] always gives [1,0,0,0]). Real panels are full of such cells; each contributed a zero to the median and the ratio inflated without bound (observed 8–211, where Hellinger distance to a centroid is bounded at 1.0, so those values were arithmetically impossible). Replaced with a Dirichlet(counts + 0.5) posterior, which represents uncertainty about *unobserved* variants: for [200,0,0,0] the truth could be [0.995,0.005,0,0]. Validated at 1.02 on diffuse and on concentrated Kv4-like nulls at variable depth, conservative (0.34) on very sparse compositions, and sensitive to real heterogeneity (1.55/3.34 mild, 3.25/8.85 strong at depth 200/2000).

*General lesson.* Compositional methods needing a pseudocount are hazardous on sparse, depth-variable data. The only reason either failure surfaced is that the noise-floor metric has a **calibration target** — pure noise must give 1.0 — so 131 was visibly impossible. The original script's R² had no such anchor, which is why its failure looked like a result. Any new metric in this project should carry a known-answer null.

**Infrastructure.** `hip_common.py` pointed `HIP_PKL`/`CTX_PKL` at paths inside an ephemeral analysis container; those pickles no longer exist anywhere. Both loaders now read the durable annotated CSVs from the collapse step (`HIP_DATA`/`CTX_DATA`). `region_compare3.py` still carries a dead scratchpad path and needs the same treatment before it can be re-run.

## Is the within-type heterogeneity biological or methodological?

Not yet settled, and the distinction matters more than anything else in this document. What the noise floor establishes is narrow but firm: **the variance is not multinomial sampling noise.** It is 4–25× that, consistently, across two chemistries. What it does *not* exclude is any source of technical variance that is depth-independent, because such a source sits in the numerator alongside biology and is invisible to this test. The √depth scaling confirms the noise floor is correctly estimated; it does **not** adjudicate the numerator.

Candidate non-biological sources, in rough order of concern:

1. **Isoform-assignment uncertainty — the main worry.** salmon's EM distributes ambiguous reads across isoforms and reports point estimates. Our analysis treats those as data, so assignment error passes straight into apparent heterogeneity. This project already has direct evidence that salmon's assignment here is imperfect: the 4-fold underestimate of Kv4.3L traced to 3′UTR model asymmetry. Short reads (50 bp PE in the Allen sets) make this worse, and variants differing only at an internal cassette exon or an alternative acceptor are exactly the hard cases.
2. **PCR amplification noise.** SMART-seq amplifies before sequencing; stochastic early-cycle amplification is multiplicative and depth-independent, so it would inflate apparent composition variance without touching the noise floor. Both v1 and v4 chemistries amplify, so cross-chemistry consistency does not rule it out.
3. **Reverse-transcription efficiency and length bias** differing between variants, varying cell to cell with RNA quality.
4. **Dissociation-induced transcriptional change.** Both protocols dissociate tissue, which induces stress and immediate-early programs and can shift splicing.
5. **Ambient RNA and index hopping**, lower in plate-based than droplet data but non-zero.

And if biological, it still splits two ways with very different consequences: **stable cell-to-cell differences** (developmental, epigenetic, allelic) versus **dynamic state** (activity, circadian, cell cycle). The first means a cell's A-current phenotype is a durable property; the second means it fluctuates in the same cell over hours.

Four tests, in order of decisiveness per unit effort:

- **salmon inferential replicates.** Re-quantify a few hundred cells with `--numBootstraps 100` (or Gibbs samples) and compare the bootstrap variance of each cell's composition against the within-subclass spread. This directly measures assignment uncertainty and would settle concern 1 outright. Cheapest decisive experiment available; the quant directories and index are already in place.
- **A negative-control gene family.** Run the identical noise-floor analysis on isoform families of comparable expression with no expected biological heterogeneity — housekeeping genes, or a structural family. If they return ratios of 10–20 as well, the variance is technical and generic. If Kv4-complex ratios are specifically high, that is strong evidence for biology. Uses existing data; a few hours.
- **Is the variance structured?** Random technical noise is uncorrelated with anything. Test whether a cell's Kv4 composition correlates with its KChIP and DPP composition, with activity/IEG state, or with dissection position. Coherent structure across the three panels within a cell would be hard to explain technically.
- **Patch-seq.** The decisive biological test, below.

Until at least the first two are done, the honest statement is: *within-type variance exceeds sampling noise by a large margin and is consistent across datasets; how much of it is isoform-assignment error is unmeasured.*

## If it is biological: predicted consequences for recorded A-currents in pyramidal neurons

The Kv4 complex is a three-component machine, and the auxiliary subunits — not the pore — set most of the biophysics. So composition variance maps onto specific, measurable current parameters rather than onto current amplitude alone.

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
- **Patch-seq is the decisive experiment**: record A-current properties, then sequence the same cell, and ask whether composition predicts kinetics. Public mouse cortical patch-seq datasets exist (Gouwens 2020, Scala 2020) and could be reanalysed before committing to new experiments — with the large caveat that patch-seq transcriptomes are shallow, and this document's own depth analysis says isoform-level resolution needs far more than they provide. That may make reanalysis infeasible and new, deeply-sequenced patch-seq the only route.

**Two caveats that could break the whole chain.** mRNA composition is not protein composition: translation efficiency, trafficking, assembly stoichiometry and protein half-life all intervene, and the lab's own evidence for direct FMRP–Kv4 protein interaction (Zhan 2024) is a reminder that regulation sits downstream of transcript abundance. And the complex assembles somatically (companion document, Oct 30), so a cell's surface channel population integrates its transcript composition over the mRNA and protein lifetime — which would *smooth* transcript-level variance and predict less current variance than the RNA suggests. If recorded A-current parameters turn out to be *tighter* than composition variance predicts, that difference is itself informative: it would imply a post-transcriptional buffering mechanism setting A-current to a target.

**Relevance to V404M.** If auxiliary-subunit composition varies this much between neighbouring pyramidal cells, the functional impact of V404M should vary with it — the mutation acts on a channel whose properties are set by partners the cell chooses. That predicts a *distribution* of severity across cells within a region rather than a uniform deficit, and offers a mechanism for why a single V404M copy is worse than one or two nulls: a null removes a channel, whereas V404M inserts a defective subunit into every complex the cell assembles, whatever its composition.

## Concerns and open questions

1. The cortex composition confound invalidates all cortical gene-level Kv4 calls until a neuronal-fraction covariate is added. Hippocampus is unaffected.
2. The XBP1 and ΔFosB junction readouts are too shallow to support ratio estimates and should not be quoted as quantitative.
3. Whether the within-type heterogeneity survives correction for isoform-assignment uncertainty is unknown and is the single most important open question here.
4. The OPC Kv4.2 observation rests on 22 cells and lacks a contamination check.
5. ΔFosB cannot be assayed properly by bulk junction counting; the IHC remains the better evidence.
6. Bulk RNA cannot exclude a small, severely reactive astrocyte subpopulation.
7. `region_compare3.py` and the three-way comparison it produced are not currently reproducible (dead paths, and the CLR metric issue applies to its distance calculations).

## Recommended additional analyses

| Priority | Analysis | Why |
| --- | --- | --- |
| 1 | salmon bootstraps on a cell subset; negative-control gene family | Decides whether the heterogeneity result stands at all |
| 2 | Neuronal-fraction covariate / deconvolution for cortex and hippocampus | Makes the V404M gene-level numbers interpretable |
| 3 | Retrospective variance analysis of existing pyramidal A-current recordings | Tests the central prediction with data already in hand |
| 4 | Neuronal-specific conditional V404M | Separates neuronal from OPC contributions *and* tests whether astrogliosis is secondary — one experiment, two questions |
| 5 | Developmental time course of the astrogliosis | Cause vs consequence |
| 6 | Neuronal-marker contamination check on Tasic 2016 OPC libraries | Cheap; closes or opens the OPC thread |
| 7 | Pdgfra⁺/Olig2⁺ counts in V404M cortex | Tests the OPC prediction directly |
| 8 | Astrocyte-resolved (sorted or snRNA) sequencing with per-cell Kcnj10 | Excludes a severely reactive subpopulation |
| 9 | Rebuild `channel_phenotype_groups.py` with the Hellinger metric and the noise floor as effect-size gate | Only worth doing after item 1 |

### Going forward

The immediate decision point is item 1. If isoform-assignment uncertainty accounts for most of the within-type spread, the heterogeneity result retracts and the by-type tables stand as they are. If it does not, then Kv4-complex composition is a per-cell property rather than a per-type one, and that changes how the survey's central tables should be read and what the ephys should be looking for — at which point item 3 becomes the cheapest way to test it against data the lab already owns.

For V404M, the biology now points away from proteotoxicity and toward circuit dysfunction: both alleles expressed, no stress response, no change in Kcnd2 processing, astrogliosis confined to epileptogenic regions and not tracking mutant dose, homeostatic astrocyte function intact. The conditional line is the experiment that would most efficiently advance this, since it addresses the neuronal-vs-OPC and the cause-vs-consequence questions at once.

## Data, code and provenance

Repository `ppfaff/Kv4-isoform-survey`. V404M pipeline: `pipeline/build_v404m_runlist.py`, `find_v404m_codon.py`, `set_v404m_haplotypes.py`, `v404m_allele_extract.py`, `run_v404m_salmon.sh`, `run_v404m_star.sh`, `collapse_v404m.py`; analysis `scripts/v404m_panels.py`, `v404m_stress.py`, `v404m_kv4.py`, `v404m_passenger.py`. Heterogeneity: `scripts/channel_noise_floor.py` (Dirichlet floor), `channel_phenotype_groups.py` (grouping, pending rebuild), `hip_common.py`. Dataset provenance for all public data in `lab_db/DATASETS.md`. Output tables in `C:\Users\paulp\kv4_v404m_tables` and `tables/` in the repository.

Companion documents: *Kv4 Complex mRNAs Are Somatic* (Oct 30) and the isoform survey conclusions (Oct 29).
