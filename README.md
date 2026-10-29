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
