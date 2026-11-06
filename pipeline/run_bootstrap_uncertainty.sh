#!/usr/bin/env bash
# Re-quantify a subset of Tasic 2016 cells with salmon inferential replicates, to measure how much of the
# apparent per-cell variant heterogeneity is isoform-ASSIGNMENT uncertainty rather than biology.
#
# Why this is needed: every per-cell analysis in this project treats salmon's point estimates as data. The EM
# apportions reads that are ambiguous between isoforms, and that apportionment has uncertainty which never
# appears in quant.sf. The Dirichlet noise floor (channel_noise_floor.py) models COUNTING noise given the
# assignment, so it cannot see this at all. The DPP6 N-terminal variants differ only at the 5' end, which is
# exactly where assignment is least certain, and they are the variants that came back bimodal.
#
# Bootstraps give the missing number directly: for each cell, the spread of a variant's fraction ACROSS
# bootstrap replicates is its assignment uncertainty. Compare that against the spread ACROSS CELLS within a
# type (variant_distribution_shape.py, column excess_sd):
#     across-cell spread >> across-bootstrap spread  -> real heterogeneity, the ephys prediction stands
#     across-cell spread ~  across-bootstrap spread  -> the heterogeneity is assignment noise, retract it
#
# Designed to run on a spare machine (32 GB RAM is ample; salmon needs a few GB for this index) and to leave
# almost nothing behind: each cell's FASTQ and its ~50 MB of raw bootstraps are deleted as soon as the
# Kv4-complex rows have been extracted, so peak disk is a few GB regardless of how many cells are run.
#
# Prerequisites:
#   micromamba create -y -n kv4boot -c conda-forge -c bioconda salmon=1.10.3 python=3.12 numpy pandas
#   micromamba activate kv4boot
#   the salmon transcriptome index (ref/salmon_idx_tx), copied from Lambda or rebuilt from GENCODE vM36
#
# Usage:
#   CELLS=cells_this_machine.txt IDX=ref/salmon_idx_tx bash run_bootstrap_uncertainty.sh
#   CELLS is one cell name per line, matching the sample_title in the Tasic 2016 ENA records
#   (e.g. the names in tables/kv4_isoform_per_cell_annotated.csv).
set -uo pipefail
CELLS=${CELLS:?set CELLS to a file of cell names, one per line}
IDX=${IDX:-ref/salmon_idx_tx}
NBOOT=${NBOOT:-100}
THREADS=${THREADS:-6}          # 4-core/8-thread i7: 6 leaves the machine usable
PROJECT=${PROJECT:-PRJNA291624}
HOST=$(hostname | tr '[:lower:]' '[:upper:]')
WORK=${WORK:-boot_work}
OUTDIR=${OUTDIR:-tables}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$WORK" "$OUTDIR" logs boot_out
log(){ echo "[$(date '+%F %T')] $*"; }
command -v salmon >/dev/null || { log "salmon not on PATH"; exit 1; }
[ -s "$IDX/versionInfo.json" ] || { log "no salmon index at $IDX"; exit 1; }

# ------------------------------------------------------------------ runlist from ENA, filtered to our cells
RUNLIST=$WORK/runlist.tsv
if [ ! -s "$RUNLIST" ]; then
  log "fetching ENA filereport for $PROJECT"
  curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${PROJECT}&result=read_run&fields=run_accession,sample_title,fastq_ftp,read_count&format=tsv&download=true" \
    > "$WORK/ena_all.tsv" || { log "ENA fetch failed"; exit 1; }
  # keep rows whose sample_title matches a requested cell name exactly
  python3 - "$WORK/ena_all.tsv" "$CELLS" > "$RUNLIST" <<'PY'
import sys, csv
want = {l.strip() for l in open(sys.argv[2]) if l.strip()}
r = csv.DictReader(open(sys.argv[1]), delimiter='\t')
w = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
w.writerow(['cell', 'run', 'fastq', 'reads'])
seen = set()
for row in r:
    t = (row.get('sample_title') or '').strip()
    if t in want and t not in seen:
        seen.add(t)
        w.writerow([t, row['run_accession'], row['fastq_ftp'], row.get('read_count', '')])
miss = want - seen
if miss:
    sys.stderr.write(f'WARNING: {len(miss)} cells not matched in ENA, e.g. {sorted(miss)[:5]}\n')
PY
fi
N=$(tail -n +2 "$RUNLIST" | wc -l)
log "$N cells to process on $HOST (requested $(grep -c . "$CELLS"))"
[ "$N" -gt 0 ] || { log "no cells matched -- check that CELLS names match ENA sample_title"; exit 1; }

# ------------------------------------------------------------------ per cell: download, quant, extract, bin
while IFS=$'\t' read -r cell run ftp reads; do
  [ -n "$cell" ] || continue
  outfile="boot_out/${cell// /_}.tsv"
  [ -s "$outfile" ] && continue                      # resumable
  d=$WORK/$run; mkdir -p "$d"
  # ENA gives one or more semicolon-separated FTP paths; single-end for this project
  ok=1
  : > "$d/files.txt"
  IFS=';' read -ra URLS <<< "$ftp"
  for u in "${URLS[@]}"; do
    f=$d/$(basename "$u")
    if [ ! -s "$f" ]; then
      curl -fsSL --retry 3 --retry-delay 5 "ftp://$u" -o "$f" || { ok=0; break; }
    fi
    echo "$f" >> "$d/files.txt"
  done
  if [ "$ok" -ne 1 ]; then
    log "$cell ($run): download failed"; echo "$cell download-failed" >> logs/boot_failures.log
    rm -rf "$d"; continue
  fi
  mapfile -t FQ < "$d/files.txt"
  salmon quant -i "$IDX" -l A -r "${FQ[@]}" -p "$THREADS" --numBootstraps "$NBOOT" \
      --seqBias --gcBias -o "$d/quant" -q > "logs/$run.salmon.log" 2>&1
  if [ -s "$d/quant/quant.sf" ] && [ -s "$d/quant/aux_info/bootstrap/bootstraps.gz" ]; then
    python3 -I "$HERE/extract_bootstraps.py" "$cell" "$d/quant" > "$outfile.part" 2>> "logs/$run.extract.log" \
      && mv "$outfile.part" "$outfile" \
      || { log "$cell: extraction failed"; echo "$cell extract-failed" >> logs/boot_failures.log; }
    log "$cell done ($(wc -l < "$outfile" 2>/dev/null || echo 0) rows)"
  else
    log "$cell ($run): quant or bootstraps missing"; echo "$cell quant-failed" >> logs/boot_failures.log
    cp "$d/quant/logs/salmon_quant.log" "logs/$run.salmon_quant.log" 2>/dev/null
  fi
  rm -rf "$d"                                        # FASTQ + raw bootstraps go immediately
done < <(tail -n +2 "$RUNLIST")

# ------------------------------------------------------------------ collect
out=$OUTDIR/bootstrap_uncertainty_${HOST}.tsv
{ head -1 "$(ls boot_out/*.tsv | head -1)"; for f in boot_out/*.tsv; do tail -n +2 "$f"; done; } > "$out"
log "wrote $out  ($(( $(wc -l < "$out") - 1 )) rows from $(ls boot_out/*.tsv | wc -l) cells)"
log "failures: $(cat logs/boot_failures.log 2>/dev/null | wc -l)"
log "copy $out back to Lambda and run: python3 scripts/bootstrap_vs_biological.py"
