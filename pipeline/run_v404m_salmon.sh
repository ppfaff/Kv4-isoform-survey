#!/usr/bin/env bash
# Salmon isoform quantification of the Kv4.2-V404M vs wild-type bulk RNA-seq (48 libraries, 4 tissues).
#
# The FASTQs live on a network drive (T:) mounted into WSL, where sustained read throughput is unpredictable, so a
# background stager copies each sample's R1/R2 to local disk a couple of samples ahead while salmon works on the
# current one. Everything is resumable: a sample with quant/<sample>/quant.sf is skipped, so the script can be
# re-run after an interruption. Staged FASTQs are deleted as soon as their quantification succeeds.
#
# Prerequisites (WSL on Lambda, kv4 env):
#   sudo mkdir -p /mnt/t && sudo mount -t drvfs 'T:' /mnt/t        # network drives are not auto-mounted
#   conda activate kv4   (salmon 1.10.3 must be on PATH)
#   python3 build_v404m_runlist.py "<csv>" "<fastq dir>"           # -> tables/v404m_samples.tsv
#
# Usage:  N_Q=2 T_Q=12 nohup nice -n 19 bash run_v404m_salmon.sh > run_v404m_salmon.log 2>&1 &
set -uo pipefail
N_Q=${N_Q:-2}              # concurrent salmon jobs
T_Q=${T_Q:-12}             # threads per job
BUFFER=${BUFFER:-2}        # samples staged ahead (each ~10 GB of FASTQ)
SAMPLES=${SAMPLES:-tables/v404m_samples.tsv}
REF_DIR=${REF_DIR:-$HOME/kv4_cortex/ref}
STAGE_DIR=${STAGE_DIR:-fastq_local}
ROOT=$(pwd)
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$STAGE_DIR" quant logs tables; rm -f "$STAGE_DIR/.all_staged"
log(){ echo "[$(date '+%F %T')] $*"; }
command -v salmon >/dev/null || { log "salmon not on PATH — conda activate kv4"; exit 1; }
[ -s "$SAMPLES" ] || { log "missing $SAMPLES — run build_v404m_runlist.py first"; exit 1; }

# ---------------------------------------------------------------- index (shared with every earlier run)
if [ -s "$REF_DIR/salmon_idx_tx/versionInfo.json" ]; then
  mkdir -p ref; [ -e ref/salmon_idx_tx ] || ln -s "$REF_DIR/salmon_idx_tx" ref/salmon_idx_tx
  log "using existing index $REF_DIR/salmon_idx_tx"
else
  log "no salmon index at $REF_DIR/salmon_idx_tx — point REF_DIR at the folder holding salmon_idx_tx"; exit 1
fi
IDX=ref/salmon_idx_tx
TOTAL=$(tail -n +2 "$SAMPLES" | wc -l); log "$TOTAL samples in $SAMPLES"

# ---------------------------------------------------------------- stager (background)
stage_one() {   # $1 sample  $2 r1  $3 r2
  s=$1; r1=$2; r2=$3
  [ -s "quant/$s/quant.sf" ] && return 0
  [ -f "$STAGE_DIR/$s/.ready" ] && return 0
  while [ "$(find "$STAGE_DIR" -name .ready 2>/dev/null | wc -l)" -ge "$BUFFER" ]; do sleep 20; done
  rm -rf "${STAGE_DIR:?}/$s"; mkdir -p "$STAGE_DIR/$s"
  if cp "$r1" "$STAGE_DIR/$s/R1.fastq.gz" && cp "$r2" "$STAGE_DIR/$s/R2.fastq.gz"; then
    touch "$STAGE_DIR/$s/.ready"
  else
    rm -rf "${STAGE_DIR:?}/$s"; echo "$s stage-failed" >> logs/v404m_failures.log
  fi
}
(
  tail -n +2 "$SAMPLES" | while IFS=$'\t' read -r s tissue genotype sex r1 r2 rest; do
    [ -n "$s" ] || continue
    stage_one "$s" "$r1" "$r2"
  done
  touch "$STAGE_DIR/.all_staged"
) &
log "stager started (pid $!, buffer $BUFFER samples)"

# ---------------------------------------------------------------- quantifier loop
quant_one() {
  s=$1
  salmon quant -i "$IDX" -l A -1 "$STAGE_DIR/$s/R1.fastq.gz" -2 "$STAGE_DIR/$s/R2.fastq.gz" \
    -p "$T_Q" --seqBias --gcBias -o "quant/$s" -q > "logs/$s.salmon.log" 2>&1
  if [ -s "quant/$s/quant.sf" ]; then
    # count the V404M allele from the staged FASTQs while they are still local: this saves a second full pass
    # over the source drive later (set ALLELE=0 to skip, e.g. before tables/v404m_codon.tsv exists)
    if [ "${ALLELE:-1}" = 1 ] && [ -s tables/v404m_codon.tsv ] && [ ! -s "allele/$s.tsv" ]; then
      mkdir -p allele
      L=$(awk -F'\t' '$1=="left_anchor"{print $2}' tables/v404m_codon.tsv)
      R=$(awk -F'\t' '$1=="right_anchor"{print $2}' tables/v404m_codon.tsv)
      LRC=$(awk -F'\t' '$1=="left_anchor_rc"{print $2}' tables/v404m_codon.tsv)
      RRC=$(awk -F'\t' '$1=="right_anchor_rc"{print $2}' tables/v404m_codon.tsv)
      zcat "$STAGE_DIR/$s/R1.fastq.gz" "$STAGE_DIR/$s/R2.fastq.gz" 2>/dev/null \
        | LC_ALL=C grep -F -e "$L" -e "$R" -e "$LRC" -e "$RRC" \
        | python3 -I "$HERE/v404m_allele_extract.py" "$s" tables/v404m_codon.tsv > "allele/$s.tsv.part" 2>> "logs/$s.allele.log"
      mv "allele/$s.tsv.part" "allele/$s.tsv"
    fi
    lt=$(grep -o '"library_types": \[[^]]*\]' "quant/$s/aux_info/meta_info.json" | tr -d '\n')
    pm=$(grep -o '"percent_mapped": [0-9.]*' "quant/$s/aux_info/meta_info.json")
    np=$(grep -o '"num_processed": [0-9]*' "quant/$s/aux_info/meta_info.json")
    echo "$s done $np $pm $lt" >> logs/v404m_progress.log
    rm -rf "${STAGE_DIR:?}/$s"
  else
    cp "quant/$s/logs/salmon_quant.log" "logs/$s.salmon_quant.log" 2>/dev/null
    rm -rf "quant/$s" "${STAGE_DIR:?}/$s"; echo "$s quant-failed $(date '+%F %T')" >> logs/v404m_failures.log
  fi
}
export -f quant_one; export IDX T_Q STAGE_DIR HERE ALLELE

last=0
while true; do
  running=$(jobs -rp | wc -l)
  if [ "$running" -lt "$N_Q" ]; then
    next=$(for d in "$STAGE_DIR"/*/; do d=${d%/}; [ -f "$d/.ready" ] && { basename "$d"; break; }; done)
    if [ -n "$next" ]; then
      mv "$STAGE_DIR/$next/.ready" "$STAGE_DIR/$next/.claimed"
      quant_one "$next" &
      continue
    fi
  fi
  if [ "$running" -eq 0 ] && [ -f "$STAGE_DIR/.all_staged" ] && [ -z "$(find "$STAGE_DIR" -name .ready -print -quit 2>/dev/null)" ]; then break; fi
  now=$(date +%s)
  if [ $((now - last)) -ge 300 ]; then
    log "quantified $(ls quant 2>/dev/null | wc -l) / $TOTAL  (running $running, staged $(find "$STAGE_DIR" -name .ready 2>/dev/null | wc -l), failures $(wc -l < logs/v404m_failures.log 2>/dev/null || echo 0))"
    last=$now
  fi
  sleep 15
done
wait
log "finished: $(ls quant | wc -l) / $TOTAL quantified; failures $(wc -l < logs/v404m_failures.log 2>/dev/null || echo 0)"
log "library types seen: $(awk '{for(i=1;i<=NF;i++) if($i ~ /^\"I/) print $i}' logs/v404m_progress.log 2>/dev/null | sort | uniq -c | tr '\n' ' ')"
if ls allele/*.tsv >/dev/null 2>&1; then
  { printf 'sample\tanchor\torientation\ttriplet\tcount\n'; cat allele/*.tsv; } > tables/v404m_allele_counts.tsv
  log "allele counts for $(ls allele/*.tsv | wc -l) samples -> tables/v404m_allele_counts.tsv"
fi
log "next: python3 collapse_v404m.py"
