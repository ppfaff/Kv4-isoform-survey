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
