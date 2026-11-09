Kv4 auxiliary-subunit isoform survey — analysis bundle (2026-10-29)

scripts/   analysis scripts (hip_common.py holds the variant map, colours and dataset loaders;
           variant_scheme.py / percell_family_genes.py / acurrent_index.py take 'hip' or 'ctx';
           region_compare3.py = three-way Tasic2016 / Allen VISp / Allen HPF comparison;
           dg_kcnd2_figure.py = DG exon vs intron; hip_lamp5_and_index.py = Lamp5 + complex index)
pipeline/  run_hip_isoforms.sh (NeMO SSv4 paired-end), run_kv4_isoforms_v5.sh (SRA/ENA single-end),
           collapse_hip.py / collapse_ctx.py (cross-machine collapse + merge), split_runlist.py,
           build_visp_runlist.py, extract_markers.py (Allen HDF5 gene extraction)
figures/   all figures delivered during the project
tables/    per-subclass / per-type tables, comparison tests, Allen marker extractions
Input per-cell tables (not included; on Lambda): ~/allen_meta/tables (HPF), ~/kv4_cortex/tables (Tasic 2016),
~/allen_visp/tables (Allen VISp). Paths inside scripts point at /home/claude/kv4/...; edit HIP_PKL / CTX_PKL
in hip_common.py and the VISp path in region_compare3.py to re-run.
Write-up: https://claude.ai/artifact/BdsYCYmxtyRcaqM968Bj9Z

Dendrite vs soma follow-up (2026-10-29) — Farris 2019 (GSE116342 / PRJNA478217, laminar LCM: cell-body vs
dendritic layers, CA1/CA2/CA3/DG x 3 animals) and Middleton 2019 (GSE115480 / PRJNA475141 / SRP150011,
16 cultured neurons, paired soma + dendrite). Run on Lambda, in a fresh work dir (e.g. ~/kv4_dendrite):
  conda activate kv4; git clone/pull this repo; cp pipeline/{build_dendrite_runlist.py,run_dendrite_isoforms.sh,collapse_dendrite.py} .
  python3 build_dendrite_runlist.py                      # tables/dendrite_samples.tsv (+ runlist); fix any '?' rows
  KEEP_FASTQ=1 N_Q=3 T_Q=8 nohup nice -n 19 bash run_dendrite_isoforms.sh > run_dendrite.log 2>&1 &
  python3 collapse_dendrite.py                           # tables/dendrite_isoform_per_sample.csv, dendrite_gene_tpm.csv.gz
  micromamba install -c bioconda -c conda-forge star subread samtools   # once
  nohup nice -n 19 bash run_dendrite_star.sh > run_star.log 2>&1 &      # Farris only: exon/intron + Kcnd2 per-intron depth
  python3 <repo>/scripts/dendrite_compartment.py tables <repo>         # figures/dendrite_*.png, tables/dendrite_*.csv
pipeline/make_intron_saf.py builds the merged-exon / intron SAF+BED for the 9 genes from the vM36 GTF
(called by run_dendrite_star.sh). KEEP_FASTQ=1 keeps the FASTQs so the STAR step does not re-download (~150 GB).
Junction check (scripts/dendrite_junctions.py): on Lambda, extract Kcnd3/Kcnip3 junctions from the STAR output with
  for f in star/*/SJ.out.tab; do r=$(basename $(dirname $f)); awk -v r=$r 'BEGIN{OFS="\t"} ($1=="chr3" && $2>105359000 && $3<105582000) || ($1=="chr2" && $2>127298000 && $3<127365000) {print r,$0}' $f; done > tables/sj_kcnd3_kcnip3.tsv
Result: Kv4.3L exon inclusion 0.17 (cell-body layers) -> 0.58 (dendritic layers), p=0.001; Kcnip3 cys-rich acceptor 0.19 -> 0.04.
Read-level Kcnd2 in DG granule cells (why are 18% of cells exonic-null but intronic-positive?):
  cd ~/kv4_dendrite && cp ~/Kv4-isoform-survey/pipeline/{build_dg_runlist.py,run_dg_star.sh} .
  python3 build_dg_runlist.py ~/allen_meta/hip_filtered_runlist.tsv ~/Kv4-isoform-survey/tables/dg_cells_markers.csv 150
  T=16 nohup nice -n 19 bash run_dg_star.sh > run_dg_star.log 2>&1 &      # re-uses ref/star_idx, genome held in shared memory
  python3 ~/Kv4-isoform-survey/scripts/dg_readlevel.py tables <repo>        # intron position / detained intron / residual junctions

V404M vs wild-type bulk RNA-seq (3-month littermates; cortex, hippocampus, thalamus, cerebellum; 48 paired-end
libraries, ~45-60 M pairs each). Sample sheet: tables/v404m_sample_info.csv. On Lambda, in ~/kv4_v404m:
  sudo mkdir -p /mnt/t && sudo mount -t drvfs 'T:' /mnt/t        # network drives are not auto-mounted in WSL
  conda activate kv4
  cp ~/Kv4-isoform-survey/pipeline/{build_v404m_runlist.py,run_v404m_salmon.sh,run_v404m_star.sh,find_v404m_codon.py,v404m_allele_extract.py,count_v404m_allele.sh,collapse_v404m.py,make_intron_saf.py} .
  mkdir -p ref tables && ln -sf ~/kv4_dendrite/ref/star_idx ref/ \
    && ln -sf ~/kv4_dendrite/ref/GRCm39.primary_assembly.genome.fa ref/ \
    && ln -sf ~/kv4_dendrite/ref/gencode.vM36.primary_assembly.annotation.gtf.gz ref/ \
    && ln -sf ~/kv4_cortex/ref/gencode.vM36.transcripts.fa.gz ref/
  cp ~/Kv4-isoform-survey/tables/v404m_sample_info.csv tables/
  python3 build_v404m_runlist.py tables/v404m_sample_info.csv "/mnt/t/Laurence David/RNA_Seq_VM_3mo"
  python3 find_v404m_codon.py                                   # -> tables/v404m_codon.tsv (V404M anchors)
  N_Q=2 T_Q=12 nohup nice -n 19 bash run_v404m_salmon.sh > run_v404m_salmon.log 2>&1 &   # ~5 h; also counts the allele
  python3 collapse_v404m.py
  N_A=2 T=12 nohup nice -n 19 bash run_v404m_star.sh > run_v404m_star.log 2>&1 &         # ~11 h
Analysis (in the container, on the small tables): scripts/v404m_stress.py (ISR/UPR/apoptosis/gliosis panels, DESeq2
per tissue) and scripts/v404m_kv4.py (allele ratio, Kv4 variants, Kcnd2 exon:intron, XBP1/dFosB/Kv4.3L junctions).
Each sample is staged to local disk before processing, so the source drive only needs ~30 MB/s to keep up.

Thalamus extension (2026-11-09) — see thalamus/README.md. Tier 1 = ABC Atlas 10x thalamus dissections +
THALMANAC MERFISH nucleus map (gene level; run on a lab machine), tier 2 = ThalamoSeq / TRN full-length
isoforms (next), tier 3 = scripts/v404m_thalamus_view.py on the existing V404M bulk tables
(tables/thalamus_bulk_kv4_summary.csv, tables/thalamus_bulk_kchip_classes.csv, figures/thalamus_bulk_kv4_complex.png).
