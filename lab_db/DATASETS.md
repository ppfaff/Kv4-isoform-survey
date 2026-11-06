# Public datasets and reference files used in the Kv4 isoform survey

Provenance for every dataset downloaded and quantified in this project, **excluding** the
two compartment datasets (Farris 2019 GSE116342 and Middleton 2019 GSE115480), which are
documented separately.

All single-cell data below is **plate-based SMART-seq**, which is why it can be used for
isoform work at all: full-length cDNA, paired-end or long single-end reads across the whole
transcript, so splice junctions and intron retention are directly observable. Droplet (10x)
data from the same Allen programmes is 3′-end-biased and cannot answer any of the questions
in this project; it was deliberately not used.

---

## 1. Tasic et al. 2016 — adult mouse primary visual cortex (pilot dataset)

| Field | Value |
|---|---|
| Accession | **GSE71585** (GEO) / **PRJNA291624** (BioProject) |
| Download route | ENA filereport API → FASTQ FTP |
| Citation | Tasic B, Menon V, Nguyen TN, et al. *Adult mouse cortical cell taxonomy revealed by single cell transcriptomics.* Nat Neurosci 19:335–346, 2016. PMID 26727548 |
| Region | Primary visual cortex (VISp/V1), Cre-driver-labelled and dissociated |
| Chemistry | SMART-Seq **v1**, **100 bp single-end** |
| Depth | ~8.7 M reads/cell (median) |
| Cells quantified | 1,677 of 1,679 annotated cells; 49 transcriptomic types |
| Role | **Pilot dataset.** Deepest per-cell coverage of anything used here, so it set the expectation for how much isoform-level signal is recoverable per cell, and established the Kv4-family variant panel before it was applied to the hippocampal data. |

Exact query used to build the runlist (`pipeline/run_kv4_isoforms_v5.sh`):

```
https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA291624&result=read_run
  &fields=run_accession,sample_accession,sample_alias,sample_title,library_layout,
          fastq_ftp,fastq_md5,read_count&format=tsv&download=true
```

Caveat: single-end 100 bp. Junction detection is fine, but there is no mate-pair constraint,
so isoform assignment for the long 3′UTR variants is weaker than in the paired-end sets.

---

## 2. Allen Institute SSv4 mouse cortex + hippocampus — hippocampal formation (main dataset)

| Field | Value |
|---|---|
| Programme | Allen Institute *Mouse Whole Cortex and Hippocampus SMART-seq*, BICCN / U19 Zeng |
| Download route | **NeMO Archive**, `https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/scell/SSv4/mouse/` — per-cell paired-end FASTQ delivered as one `.fastq.tar` per cell |
| Citation | Yao Z, van Velthoven CTJ, Nguyen TN, et al. *A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation.* Cell 184:3222–3241, 2021. PMID 34004146 |
| Regions kept | HIP, HIP-CA, HIPa, HIPma, HIPmp, HIPp, SUB-ProS (plus adjoining cortex in the parent release) |
| Chemistry | SMART-Seq **v4**, **50 bp paired-end** |
| Depth | ~0.70 M reads/cell (median) — **≈12× shallower than Tasic 2016** |
| Cells quantified | 8,131 of 8,137 matched hippocampal-formation cells (~390 GB raw); 26 subclasses, 18 with ≥15 cells |
| Role | **Main regional comparison.** All CA1/CA2/CA3/DG/SUB Kv4-variant composition, family co-occurrence, A-current index and Kv4-null-cell results come from this set. |

Pipeline: `pipeline/run_hip_isoforms.sh` (download → extract → salmon → delete, per cell),
runlist filtered to hippocampal formation only.

Notes on this set:
- NeMO returns HTTP errors with a 200-like body under load, so the downloader requires
  `curl -f` and size verification; without that, truncated tars masquerade as corrupt cells.
- The 50 bp paired-end reads are short. Junction-spanning reads need ≥8 nt overhang on each
  side, so PSI estimates are noisier per cell than in the compartment datasets, which is why
  all variant fractions here are aggregated by subclass rather than reported per cell.
- The GEO mirror of this programme was **not** used; NeMO was the download source, and the
  cell barcodes/sample names in `tables/hip_*.csv` are NeMO/Allen specimen names.

---

## 3. Allen Institute SSv4 VISp subset (chemistry control)

| Field | Value |
|---|---|
| Programme | Same Allen SSv4 release as §2, VISp arm (the cells reported in Tasic et al. 2018) |
| Download route | `https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/scell/SSv4/mouse/raw/VISp/` (directory listing parsed by `pipeline/build_visp_runlist.py`) |
| Citation | Tasic B, Yao Z, Graybuck LT, et al. *Shared and distinct transcriptomic cell types across neocortical areas.* Nature 563:72–78, 2018. PMID 30382198 |
| Region | Primary visual cortex (VISp/V1) |
| Chemistry | SMART-Seq **v4**, **50 bp paired-end** |
| Depth | ~2.5 M reads/cell |
| Cells quantified | 1,511 cells — a **stratified subset, 100 cells per subclass**, not the full release |
| Role | **Chemistry control.** Same brain region as §1 but the same chemistry and read length as §2. Any difference between §1 and §2 could be region *or* platform; §3 separates the two. Variant fractions in §3 track §1 rather than being pulled toward §2, which is what licenses the regional interpretation of the hippocampal results. |

This is the single most important control in the project and the reason the Tasic 2016 /
hippocampus comparison can be read as biology rather than as a depth artefact.

---

## 4. Allen combined exon+intron count matrix (marker genes only)

| Field | Value |
|---|---|
| File | `expression_matrix.hdf5` |
| URL | `https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_smart-seq/expression_matrix.hdf5` |
| Content | Allen's own released **exon + intron** gene-level counts for the SSv4 cortex+hippocampus release (same cells as §2/§3) |
| Used by | `pipeline/extract_markers.py` → `tables/allen_markers.csv`, `tables/allen_exon_markers.csv` |
| Role | **Cell-type marker QC only.** Not used for any Kv4 quantification. |

Why a second source for the same cells: the Allen matrix counts reads in introns as well as
exons, whereas our salmon quantification is transcriptome-only. Comparing the two gives an
independent per-cell check that a cell is the subclass its label claims (marker genes behave
the same way in both), and the exon-only vs exon+intron pair is also what first flagged how
much intronic signal these libraries carry — the observation that led to the exon:intron
ratio readout used throughout the project.

Caveat: these are Allen's counts from Allen's own alignment, not ours. They are used for
marker QC and never mixed into a Kv4 number.

---

## 5. Reference files (GENCODE vM36 / GRCm39)

Base URL: `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36`

| File | Use |
|---|---|
| `gencode.vM36.transcripts.fa.gz` | salmon selective-alignment index (`ref/salmon_idx_tx`), transcriptome-only |
| `gencode.vM36.primary_assembly.annotation.gtf.gz` | STAR index, featureCounts, exon/intron SAF construction, codon lookup |
| `GRCm39.primary_assembly.genome.fa.gz` | STAR genome, `samtools faidx` for the Kcnd2 CDS and the V404M codon |

One annotation release is used for every dataset in the project, including the V404M bulk
RNA-seq, so isoform IDs and junction coordinates are directly comparable across all of them.
All genomic coordinates quoted anywhere in this project are **GRCm39**, not mm10 — the Kcnd2
locus differs by roughly 0.1 Mb between the two assemblies, so coordinates from older
literature must be lifted before comparison.

Tool versions: salmon 1.10.3, STAR 2.7.11b, subread (featureCounts) 2.1.1, samtools 1.24,
pydeseq2 0.5.4.

---

## What is *not* in this list

- **Farris 2019 (GSE116342)** and **Middleton 2019 (GSE115480)** — the soma/neuropil and
  soma/dendrite compartment datasets. Documented in `KV4_DENDRITE_COMPARTMENT.md`.
- **V404M bulk RNA-seq** — unpublished in-house data (48 libraries, 4 tissues, littermate
  pairs), not downloaded.
- Any 10x/droplet Allen data — excluded as 3′-biased, as above.
- No miRNA-seq dataset has been downloaded yet; the cell-type-specific miRNA question raised
  earlier remains open and would need a new acquisition.

---

*Compiled from the download code in this repository (`pipeline/`), which is the authoritative
record of what was actually fetched.*

---

## Candidate sources for DIRECT mRNA half-life measurement (not yet downloaded)

The exon:intron estimate in `KV4_V404M_BULK_AND_HETEROGENEITY.md` infers Kcnd2's half-life as
2.6-13 h with the splicing commitment time as the dominant uncertainty. Two published datasets
measure decay directly and would replace that inference.

### 1. Qiu, Zhang, Xia, Gao et al., *Nat Neurosci* 29(11):2685-2698 (Nov 2026). PMID 42675226

In vivo transgenesis-free 4sU labelling in mouse brain with SLAM chemistry; droplet scNT-seq2 plus
spatial NT-seq on Visium. Reports gene-level half-lives for **12,483 genes across 15 cortical cell
types**, mean 4.6 h, range <1 h to >100 h; neurite-enriched mRNAs mean 11.9 h. Our 112:1 density
ratio for Kcnd2 puts it at 2.6-13 h, i.e. squarely on their transcriptome average.

| accession | content | size | format |
|---|---|---|---|
| GSE249405 | SuperSeries | 6.7 Gb | — |
| GSE249398 | scNT-seq2, cortical cell types | 2.0 Gb | **tar of H5AD** |
| GSE249403 | Spatial NT-seq | 4.6 Gb | tar |
| PRJNA1048969 | BioProject | — | — |

The H5AD files carry per-cell **new and total** RNA counts, so no FASTQ, no alignment and no
conversion-calling is needed. Pseudobulk new and total counts within a cell type, then
`tau = -t / ln(1 - NTR)` and `half-life = tau * ln2` at the known labelling time t. Include *Fos*
(minutes) and *Actb*/*Gapdh* (hours) as known-answer anchors before believing any Kv4-panel number.

**Caveat that decides whether this is worth it: the labelled animals are E16.5, P7 and P12** —
developing brain, not adult. Kcnd2 is developmentally regulated, so these are a weaker reference for
adult V404M recordings than the cell-type resolution suggests.

Supplementary Table S2 (per-gene synthesis, degradation and abundance across the 15 cortical types)
would answer the panel question on its own without the 2 Gb download. It is **not** in GEO and
**not** attached to the bioRxiv/PMC preprint (PMC12874000 carries the PDF only; its availability
statement still says "upon publication"), so it has to come from the *Nat Neurosci* article.

### 2. Loedige et al., *Mol Cell* (2023) — mRNA stability and m6A in neuronal localization

Bulk SLAM-seq pulse-chase in **mouse primary cortical neurons** and mESC-derived neurons, with
half-lives resolved by compartment (soma vs neurites). Bulk gives far better per-gene precision than
single-cell NTR, and the compartment split parallels the Farris soma/dendrite libraries already used
here. Deposited in ArrayExpress, **E-MTAB-9400** and **E-MTAB-12888**; per-gene values in its
Tables S1/S2, codon optimality vs half-life in S7. Reports that longer 3'UTRs are destabilising,
which is a direct testable prediction for Kcnd2's long 3'UTR.

### Access note

Neither the analysis container nor the linked Windows VM can reach ncbi.nlm.nih.gov, ebi.ac.uk or
biorxiv.org — both sit behind an egress allowlist that passes only package registries and GitHub
(the VM does reach pypi, so `h5py`/`anndata` can be installed there). The desktop browser pane
reaches all of them but each host needs a separate approval, and it has 1.9 Gb free on its session
disk, so a multi-Gb download must land in a connected Windows folder and be processed in place.
