#!/usr/bin/env bash
# Salmon isoform quantification for the dendrite-vs-soma datasets (Farris 2019 + Middleton 2019).
# Same decoupled, resumable design as run_kv4_isoforms_v5.sh: a background downloader fills a bounded
# FASTQ buffer (SRA prefetch + fasterq-dump, falling back to ENA's FASTQ mirror) while a streaming
# quantifier keeps N_Q salmon jobs busy. Bulk libraries here are 20-50 M read pairs each, so the
# defaults are fewer, bigger jobs than the single-cell runs.
#
# Set-up (once):   python3 build_dendrite_runlist.py        -> tables/dendrite_runlist.tsv (56 runs)
# Usage:           N_Q=3 T_Q=8 nohup nice -n 19 bash run_dendrite_isoforms.sh > run_dendrite.log 2>&1 &
# Re-use the existing index:  REF_DIR=$HOME/kv4_cortex/ref  (default below; it is built here if absent)
# KEEP_FASTQ=1 leaves the FASTQs in place after quantification (marked fastq/<run>/.quantified) so that
# run_dendrite_star.sh can align the Farris libraries without downloading them again (~150 GB kept).
set -uo pipefail
KEEP_FASTQ=${KEEP_FASTQ:-0}
N_DL=${N_DL:-2}          # parallel downloads (SRA is happy with 2-3)
N_Q=${N_Q:-3}            # concurrent salmon jobs
T_Q=${T_Q:-8}            # threads per salmon job
BUFFER=${BUFFER:-6}      # downloaded-but-unquantified libraries to hold (each ~5-10 GB of FASTQ)
RUNLIST=${RUNLIST:-tables/dendrite_runlist.tsv}
REF_DIR=${REF_DIR:-$HOME/kv4_cortex/ref}
ROOT=$(pwd)
mkdir -p fastq quant tables logs sra; rm -f fastq/.all_downloaded
log(){ echo "[$(date '+%F %T')] $*"; }

[ -s "$RUNLIST" ] || { log "missing $RUNLIST — run build_dendrite_runlist.py first"; exit 1; }
TOTAL=$(wc -l < "$RUNLIST"); log "$TOTAL runs in $RUNLIST"

# ---------------------------------------------------------------- 1. transcriptome index (shared with the cortex/hip runs)
GC=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36
if [ -s "$REF_DIR/salmon_idx_tx/versionInfo.json" ]; then
  mkdir -p ref
  [ -e ref/salmon_idx_tx ] || ln -s "$REF_DIR/salmon_idx_tx" ref/salmon_idx_tx
  [ -e ref/gencode.vM36.transcripts.fa.gz ] || ln -s "$REF_DIR/gencode.vM36.transcripts.fa.gz" ref/gencode.vM36.transcripts.fa.gz
  log "using existing index $REF_DIR/salmon_idx_tx"
else
  mkdir -p ref; cd ref
  [ -s gencode.vM36.transcripts.fa.gz ] || curl -sS -O $GC/gencode.vM36.transcripts.fa.gz
  if [ ! -s salmon_idx_tx/versionInfo.json ]; then
    log "building transcriptome-only index"
    salmon index -t gencode.vM36.transcripts.fa.gz -i salmon_idx_tx -k 31 --gencode -p 8 > ../logs/index_tx.log 2>&1 || { log "index failed, see logs/index_tx.log"; exit 1; }
  fi
  cd "$ROOT"
fi
IDX=ref/salmon_idx_tx

# ---------------------------------------------------------------- 2. downloader (background)
if command -v prefetch >/dev/null 2>&1 && command -v fasterq-dump >/dev/null 2>&1; then DL_MODE=sra; else DL_MODE=ena; fi
log "download mode: $DL_MODE"

download_sra() {   # $1 = run
  run=$1; got=0
  for attempt in 1 2 3; do
    timeout -k 30 5400 prefetch -q --max-size u -O sra "$run" >> logs/$run.prefetch.log 2>&1 && { got=1; break; }
    echo "attempt $attempt timed out/failed" >> logs/$run.prefetch.log; rm -rf sra/$run; sleep 60
  done
  [ $got = 1 ] \
    && fasterq-dump -q -e 4 --split-3 -t sra/$run/tmp -O fastq/$run sra/$run > logs/$run.dump.log 2>&1 \
    && ls fastq/$run/*.fastq >/dev/null 2>&1
  rc=$?
  rm -rf sra/$run
  if [ $rc -eq 0 ] && [ -s fastq/$run/${run}_1.fastq ] && [ -s fastq/$run/${run}_2.fastq ]; then rm -f fastq/$run/${run}.fastq; fi
  return $rc
}

download_ena() {   # $1 = run, $2 = fastq_ftp (one or two ';'-separated paths)
  run=$1; ftp=$2
  [ -n "$ftp" ] || return 1
  IFS=';' read -ra parts <<< "$ftp"
  for p in "${parts[@]}"; do
    f=$(basename "$p"); ok=0
    for attempt in 1 2 3 4; do
      curl -fsS --retry 2 --retry-delay 20 -o "fastq/$run/$f" "https://$p" && { ok=1; break; }
      echo "$run $f attempt $attempt failed" >> logs/retries.log; sleep $((attempt * 20))
    done
    [ $ok = 1 ] || return 1
    gzip -t "fastq/$run/$f" 2>/dev/null || { echo "$run $f gzip-corrupt" >> logs/retries.log; return 1; }
  done
  return 0
}

download_one() {
  run=$1; ftp=${2:-}
  [ -s quant/$run/quant.sf ] && return 0
  [ -f fastq/$run/.ready ] && return 0
  while [ "$(ls fastq 2>/dev/null | while read d; do [ -f fastq/$d/.ready ] && echo x; done | wc -l)" -ge "$BUFFER" ]; do sleep 30; done
  rm -rf fastq/$run; mkdir -p fastq/$run
  if [ "$DL_MODE" = sra ]; then download_sra "$run"; else download_ena "$run" "$ftp"; fi
  if [ $? -eq 0 ] && ls fastq/$run/*.f*q* >/dev/null 2>&1; then
    touch fastq/$run/.ready
  else
    rm -rf fastq/$run; echo "$run download-failed" >> logs/failures.log
  fi
}
export -f download_one download_sra download_ena; export BUFFER DL_MODE
( awk -F'\t' '{print $1"\t"$2}' "$RUNLIST" | xargs -P "$N_DL" -L 1 bash -c 'download_one "$0" "$1"'; touch fastq/.all_downloaded ) &
log "downloader started (pid $!, $N_DL parallel, mode $DL_MODE)"

# ---------------------------------------------------------------- 3. quantifier loop
quant_one() {
  run=$1
  files=( $(ls fastq/$run/*.f*q* | grep -v '\.ready\|\.claimed') )
  if [ ${#files[@]} -eq 2 ]; then
    salmon quant -i "$IDX" -l A -1 "${files[0]}" -2 "${files[1]}" -p "$T_Q" --validateMappings --seqBias --gcBias -o quant/$run -q > logs/$run.log 2>&1
  else
    salmon quant -i "$IDX" -l A -r "${files[0]}" -p "$T_Q" --validateMappings --seqBias -o quant/$run -q > logs/$run.log 2>&1
  fi
  if [ -s quant/$run/quant.sf ]; then
    if [ "${KEEP_FASTQ:-0}" = 1 ]; then touch fastq/$run/.quantified; else rm -rf fastq/$run; fi
    echo "$run done $(grep -o '"percent_mapped": [0-9.]*' quant/$run/aux_info/meta_info.json) $(grep -o '"library_types": \[[^]]*\]' quant/$run/aux_info/meta_info.json)" >> logs/progress.log
  else
    rm -rf fastq/$run quant/$run; echo "$run quant-failed" >> logs/failures.log
  fi
}
export -f quant_one; export IDX T_Q KEEP_FASTQ
last_log=$(date +%s)
while true; do
  running=$(jobs -rp | wc -l)
  if [ "$running" -lt "$N_Q" ]; then
    next=$(for d in fastq/*/; do d=${d%/}; [ -f "$d/.ready" ] && { echo "${d#fastq/}"; break; }; done)
    if [ -n "$next" ]; then
      mv "fastq/$next/.ready" "fastq/$next/.claimed"
      quant_one "$next" &
      continue
    fi
  fi
  if [ "$running" -eq 0 ] && [ -f fastq/.all_downloaded ] && [ -z "$(find fastq -name .ready -print -quit 2>/dev/null)" ]; then break; fi
  now=$(date +%s)
  if [ $((now - last_log)) -ge 300 ]; then
    log "quantified so far: $(ls quant | wc -l) / $TOTAL  (running: $running, ready: $(find fastq -name .ready 2>/dev/null | wc -l))"
    last_log=$now
  fi
  sleep 10
done
wait
log "quantification finished: $(ls quant | wc -l) runs; failures: $(wc -l < logs/failures.log 2>/dev/null || echo 0)"
log "next: python3 collapse_dendrite.py  (then scripts/dendrite_compartment.py)"
