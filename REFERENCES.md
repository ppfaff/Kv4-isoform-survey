# Datasets analysed in the Kv4 isoform survey — references and accessions

Verified against publisher / PubMed / GEO records on 2026-11-14 except where marked.

## Cortex (variant-level, full-length SMART-seq)

**Tasic 2016 — primary visual cortex, 1,679 cells, 49 transcriptomic types.**
Tasic B, Menon V, Nguyen TN, Kim TK, Jarsky T, Yao Z, Levi B, Gray LT, Sorensen SA, Dolbeare T, Bertagnolli D, Goldy J, Shapovalova N, Parry S, Lee C, Smith K, Bernard A, Madisen L, Sunkin SM, Hawrylycz M, Koch C, Zeng H. Adult mouse cortical cell taxonomy revealed by single cell transcriptomics. *Nat Neurosci* 2016;19(2):335–346. doi:10.1038/nn.4216. GEO GSE71585.
Used for: cortical N-terminal variant fractions by type (`tables/ctx_variant_*`), cortical loop populations (`tables/ctx_loop_types_*`), the interneuron proxy for TRN variants, the bulk-vs-single-cell calibration.

**Allen Institute SMART-seq v4, isocortex + hippocampal formation (Yao 2021).**
Yao Z, van Velthoven CTJ, Nguyen TN, Goldy J, Sedeno-Cortes AE, Baftizadeh F, Bertagnolli D, Casper T, Chiang M, Crichton K, Ding S-L, Fong O, Garren E, Glandon A, Gouwens NW, Gray J, Graybuck LT, Hawrylycz MJ, Hirschstein D, Kroll M, Lathia K, Lee C, Levi B, McMillen D, Mok S, Pham T, Ren Q, Rimorin C, Shapovalova N, Sulc J, Sunkin SM, Tieu M, Torkelson A, Tung H, Ward K, Dee N, Smith KA, Tasic B, Zeng H. A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation. *Cell* 2021;184(12):3222–3241.e26. doi:10.1016/j.cell.2021.04.021. Data: NeMO dat-jb2f34y; GEO GSE185862 (Mouse Whole Cortex and Hippocampus SMART-seq); portal.brain-map.org.
Used for: hippocampal (CA1, CA3, DG, SUB) and additional cortical variant fractions (`tables/hip_variant_*`), CA1/DG reference rows of the loop figure.

## Dendrite vs soma (variant-level)

**Farris 2019 — laminar LCM, CA1/CA2/CA3/DG.**
Farris S, Ward JM, Carstens KE, Samadi M, Wang Y, Dudek SM. Hippocampal subregions express distinct dendritic transcriptomes that reveal differences in mitochondrial function in CA2. *Cell Rep* 2019;29(2):522–539.e6. doi:10.1016/j.celrep.2019.08.093. GEO GSE116343.

**Middleton 2019 — single-neuron soma/dendrite.**
Middleton SA, Eberwine J, Kim J. Comprehensive catalog of dendritically localized mRNA isoforms from sub-cellular sequencing of single mouse neurons. *BMC Biol* 2019;17:5. doi:10.1186/s12915-019-0630-z. GEO GSE115480 (SRA SRP150011).

## Thalamus

**Allen Brain Cell Atlas, whole mouse brain 10x (Yao 2023) — gene level.**
Yao Z, van Velthoven CTJ, Kunst M, Zhang M, McMillen D, Lee C, Jung W, Goldy J, Abdelhak A, Aitken M, Baker K, Baker P, Barkan E, Bertagnolli D, Bhandiwad A, Bielstein C, Bishwakarma P, Campos J, Carey D, Casper T, Chakka AB, Chakrabarty R, Chavan S, Chen M, Clark M, Close J, Crichton K, Daniel S, DiValentin P, Dolbeare T, Ellingwood L, Fiabane E, Faber T, Fliss T, Gee J, Gerstenberger J, Glandon A, Gloe J, Gould J, Gray J, Guilford N, Guzman J, Hirschstein D, Ho W, Hooper M, Huang M, Hupp M, Jin K, Kroll M, Lathia K, Leon A, Li S, Long B, Madigan Z, Malloy J, Malwane J, Martin R, McCue R, McGinty R, Mei N, Melchor J, Meyerdierks E, Mollenkopf T, Moonsman S, Nguyen TN, Otto S, Pham T, Rimorin C, Ruiz A, Sanchez R, Sawyer L, Shapovalova N, Shepard N, Slaughterbeck C, Sulc J, Tieu M, Torkelson A, Tung H, Valera Cuevas N, Vance S, Wadhwani K, Ward K, Levi B, Farrell C, Young R, Staats B, Wang MM, Thompson CL, Mufti S, Pagan CM, Kruse L, Dee N, Sunkin SM, Esposito L, Hawrylycz MJ, Waters J, Ng L, Smith K, Tasic B, Zhuang X, Zeng H. A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain. *Nature* 2023;624(7991):317–332. doi:10.1038/s41586-023-06812-z. Data: Allen Brain Cell Atlas (alleninstitute.github.io/abc_atlas_access), release 20250531; 10x matrices WMB-10Xv2/v3 TH, HY, MB.
Used for: thalamic cell-type map of the Kv4 complex (`tables/th10x_*`), TRN core/shell.

**Allen Brain Cell Atlas, whole-brain MERFISH (Zhang 2023) — spatial reference used by THALMANAC.**
Zhang M, Pan X, Jung W, Halpern AR, Eichhorn SW, Lei Z, Cohen L, Smith KA, Tasic B, Yao Z, Zeng H, Zhuang X. Molecularly defined and spatially resolved cell atlas of the whole mouse brain. *Nature* 2023;624(7991):343–354. doi:10.1038/s41586-023-06808-9.
(Author list abbreviated here; verify against the article.)

**THALMANAC (thalamus subset, nucleus annotations, analysis package).**
Turner MA, Chartrand T, Summers MT, Hooper M, van Velthoven C, Waters J, de Vries S, Zeng H, Tasic B, Svoboda K, Long B. Exploring the correspondence between gene expression and thalamic nuclei using the THALMANAC resource. *bioRxiv* 2025, doi:10.1101/2025.09.30.679413 (posted 2 Oct 2025). Code: github.com/AllenNeuralDynamics/thalmanac and github.com/AllenNeuralDynamics/abc-merfish-analysis (package `abc-merfish-analysis`).
Used for: taxonomy slice and cluster→nucleus annotations (`thalamus/annotations_c2n_combined.csv`, `thalamus/taxonomy_CCN20230722_th_zi_only.csv`).

**ThalamoSeq (Phillips 2019) — relay-neuron variant level, full-length pooled cells.**
Phillips JW, Schulmann A, Hara E, Winnubst J, Liu C, Valakh V, Wang L, Shields BC, Korff W, Chandrashekar J, Lemire AL, Mensh B, Dudman JT, Nelson SB, Hantman AW. A repeated molecular architecture across thalamic pathways. *Nat Neurosci* 2019;22(11):1925–1935. doi:10.1038/s41593-019-0483-3. GEO GSE133911 (pooled, 120 samples; SRA SRP213359), GSE133912 (single cell). Browser: thalamoseq.janelia.org.
(Middle authors from memory; verify against the article.)

**TRN subnetworks (Hartley 2024) — Spp1⁺ / Ecel1⁺ definition.**
Hartley ND, Krol A, Choi S, Rome N, Levandowski K, Pasqualoni S, Jones C, Tian J, Lee S, Lee H, Kast R, Feng G, Fu Z. Distinct structural and functional connectivity of genetically segregated thalamoreticular subnetworks. *Cell Rep* 2024;43(12):115037. doi:10.1016/j.celrep.2024.115037.
(Volume/issue/article number from memory; verify.)

## In-house

**Kv4.2-V404M/+ vs +/+ bulk RNA-seq** — 3-month littermate pairs, 48 polyA libraries (cortex, hippocampus, thalamus, cerebellum; 3′-biased), STAR + salmon, GENCODE vM36. Unpublished; lab database documents 12, 18 and the Nov 2026 revisions.

## Reference genome / annotation
GENCODE vM36 (GRCm39) transcript set for all salmon quantification; variant nomenclature per Jerng & Pfaffinger (2014 review) with the lab's extensions (KChIP1e, KChIP2b, Kcnip3 Cys-rich, KChIP4d).
