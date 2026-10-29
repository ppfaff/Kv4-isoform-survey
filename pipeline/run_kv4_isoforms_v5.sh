#!/usr/bin/env bash
# Kv4 isoform quantification for GSE71585 / PRJNA291624 — v5
# Same design as v4 (decoupled downloader + salmon quantifier, resumable), with two changes so the
# remainder of a run can be handed to another machine:
#   * RUNLIST is overridable (default: tables/runlist.tsv, regenerated from ENA only if missing).
#     Build a "remaining cells" list with awk and pass RUNLIST=tables/runlist_sweep.tsv.
#   * Downloader uses NCBI prefetch+fasterq-dump when the SRA toolkit is installed, and otherwise
#     falls back to ENA's direct FASTQ mirror (column 2 of the runlist) with curl -f, size-checked.
# Usage: nohup bash run_kv4_isoforms_v5.sh > run.log 2>&1 &
set -uo pipefail
N_DL=${N_DL:-8}
N_Q=${N_Q:-2}
T_Q=${T_Q:-4}
BUFFER=${BUFFER:-12}
ROOT=$(pwd)
mkdir -p ref fastq quant tables logs sra; rm -f fastq/.all_downloaded   # stale marker from a finished run would make the quantifier exit immediately
log(){ echo "[$(date '+%F %T')] $*"; }

# ---------------------------------------------------------------- 1. run table
RUNLIST=${RUNLIST:-tables/runlist.tsv}
if [ ! -s tables/ena_runs.tsv ]; then
  curl -sS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA291624&result=read_run&fields=run_accession,sample_accession,sample_alias,sample_title,library_layout,fastq_ftp,fastq_md5,read_count&format=tsv&download=true" > tables/ena_runs.tsv
fi
if [ "$RUNLIST" = "tables/runlist.tsv" ] && [ ! -s tables/runlist.tsv ]; then
  tail -n +2 tables/ena_runs.tsv | awk -F'\t' '$6!=""{print $1"\t"$6}' > tables/runlist.tsv
fi
[ -s "$RUNLIST" ] || { log "missing $RUNLIST"; exit 1; }
TOTAL=$(wc -l < "$RUNLIST"); log "$TOTAL runs in $RUNLIST"

# ---------------------------------------------------------------- 2. transcriptome-only index
cd ref
GC=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36
[ -s gencode.vM36.transcripts.fa.gz ] || curl -sS -O $GC/gencode.vM36.transcripts.fa.gz
if [ ! -s salmon_idx_tx/versionInfo.json ]; then
  log "building transcriptome-only index"
  salmon index -t gencode.vM36.transcripts.fa.gz -i salmon_idx_tx -k 31 --gencode -p 8 > ../logs/index_tx.log 2>&1 || { log "index failed, see logs/index_tx.log"; exit 1; }
fi
cd "$ROOT"
IDX=ref/salmon_idx_tx

# ---------------------------------------------------------------- 3. downloader (background)
if command -v prefetch >/dev/null 2>&1 && command -v fasterq-dump >/dev/null 2>&1; then
  DL_MODE=sra
else
  DL_MODE=ena
fi
log "download mode: $DL_MODE"

download_sra() {   # $1 = run
  run=$1; got=0
  for attempt in 1 2 3; do
    timeout -k 30 2400 prefetch -q --max-size u -O sra "$run" >> logs/$run.prefetch.log 2>&1 && { got=1; break; }
    echo "attempt $attempt timed out/failed" >> logs/$run.prefetch.log; rm -rf sra/$run; sleep 60
  done
  [ $got = 1 ] \
    && fasterq-dump -q -e 2 --split-3 -t sra/$run/tmp -O fastq/$run sra/$run > logs/$run.dump.log 2>&1 \
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
  rm -rf fastq/$run; mkdir -p fastq/$run   # start clean: leftovers from an interrupted download would be mistaken for a mate file
  if [ "$DL_MODE" = sra ]; then download_sra "$run"; else download_ena "$run" "$ftp"; fi
  if [ $? -eq 0 ] && ls fastq/$run/*.f*q* >/dev/null 2>&1; then
    touch fastq/$run/.ready
  else
    rm -rf fastq/$run; echo "$run download-failed" >> logs/failures.log
  fi
}
export -f download_one download_sra download_ena; export BUFFER DL_MODE
( awk -F'\t' '{print $1"\t"$2}' "$RUNLIST" | xargs -P "$N_DL" -L 1 bash -c 'download_one "$0" "$1"'; touch fastq/.all_downloaded ) &
DL_PID=$!
log "downloader started (pid $DL_PID, $N_DL parallel, mode $DL_MODE)"

# ---------------------------------------------------------------- 4. quantifier loop
quant_one() {
  run=$1
  files=( $(ls fastq/$run/*.f*q* | grep -v '\.ready') )
  if [ ${#files[@]} -eq 2 ]; then
    salmon quant -i "$IDX" -l A -1 "${files[0]}" -2 "${files[1]}" -p "$T_Q" --validateMappings --seqBias --gcBias -o quant/$run -q > logs/$run.log 2>&1
  else
    salmon quant -i "$IDX" -l A -r "${files[0]}" -p "$T_Q" --validateMappings --seqBias -o quant/$run -q > logs/$run.log 2>&1
  fi
  if [ -s quant/$run/quant.sf ]; then
    rm -rf fastq/$run; echo "$run done $(grep -o '"percent_mapped": [0-9.]*' quant/$run/aux_info/meta_info.json)" >> logs/progress.log
  else
    rm -rf fastq/$run; rm -rf quant/$run; echo "$run quant-failed" >> logs/failures.log
  fi
}
export -f quant_one; export IDX T_Q
# Streaming scheduler: keep N_Q salmon jobs running at all times, starting a new one the moment
# any finishes. (The earlier batch-per-pass design waited for the slowest cell in each batch,
# which left most workers idle whenever cell depth is uneven — as it is for Tasic 2016.)
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
  sleep 5
done
wait
log "quantification finished: $(ls quant | wc -l) runs; failures: $(wc -l < logs/failures.log 2>/dev/null || echo 0)"
log "all done — run collapse_ctx.py to build tables (merges across machines)"
