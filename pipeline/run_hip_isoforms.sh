#!/usr/bin/env bash
# Kv4 isoform quantification for Allen Institute SSv4 mouse cortex+hippocampus SMART-seq
# (Yao et al. 2021 / NeMO BICCN), filtered to hippocampal-formation cells only
# (HIP, HIP-CA, HIPa, HIPma, HIPmp, HIPp, SUB-ProS = 8,137 matched cells, ~390GB raw).
#
# Reuses the same transcriptome-only salmon index / download-quantify-delete design as
# run_kv4_isoforms_v4.sh, but downloads per-cell paired-end FASTQ as a .fastq.tar from
# NeMO (not SRA), extracts the R1/R2 pair, quantifies, then deletes.
#
# Requires: a runlist tsv (from join_hip.py, or a split_runlist.py chunk) in the working directory.
# Usage: cd ~/allen_meta && RUNLIST=hip_runlist_<machine>.tsv nohup bash run_hip_isoforms.sh > run.log 2>&1 &
# (RUNLIST defaults to hip_filtered_runlist.tsv — the full, unsplit list — if not set.)
set -uo pipefail
N_DL=${N_DL:-8}
N_Q=${N_Q:-2}
T_Q=${T_Q:-4}
BUFFER=${BUFFER:-16}
ROOT=$(pwd)
mkdir -p ref fastq quant tables logs; rm -f fastq/.all_downloaded   # stale marker from a finished run would make the quantifier exit immediately
log(){ echo "[$(date '+%F %T')] $*"; }

RUNLIST=${RUNLIST:-hip_filtered_runlist.tsv}
[ -s "$RUNLIST" ] || { log "missing $RUNLIST — run join_hip.py (or split_runlist.py) first"; exit 1; }
TOTAL=$(($(wc -l < "$RUNLIST") - 1))
log "$TOTAL hippocampal cells to process"

# ---------------------------------------------------------------- 1. transcriptome-only index
cd ref
GC=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36
[ -s gencode.vM36.transcripts.fa.gz ] || curl -sS -O $GC/gencode.vM36.transcripts.fa.gz
if [ ! -s salmon_idx_tx/versionInfo.json ]; then
  log "building transcriptome-only index"
  salmon index -t gencode.vM36.transcripts.fa.gz -i salmon_idx_tx -k 31 --gencode -p 8 > ../logs/index_tx.log 2>&1 \
    || { log "index failed, see logs/index_tx.log"; exit 1; }
fi
cd "$ROOT"
IDX=ref/salmon_idx_tx

# ---------------------------------------------------------------- 2. downloader (background)
# verifies downloaded size against the runlist's expected size (col 2) before trusting a
# transfer, and retries with backoff. Critically, curl must be told -f/--fail: without it,
# curl exits 0 on an HTTP error response (e.g. NeMO returning "503 Service Temporarily
# Unavailable" under heavy combined load) and happily saves the HTML error page as if it
# were the real file — that masquerades as corruption (tar-extract-failed/missing-R1R2)
# when it's actually just an unhandled server error that should have been retried.
#
# Circuit breaker: all download workers share a consecutive-failure counter. Once it reaches
# BREAK_AFTER, workers stop issuing requests and one of them probes a known-good URL every
# PAUSE_SECS; downloads resume only after the probe succeeds. Without this, an outage makes the
# script churn through the entire runlist in a few hours, logging thousands of failures while
# keeping the request rate pinned — the worst possible behaviour toward a struggling server.
BREAK_AFTER=${BREAK_AFTER:-10}
PAUSE_SECS=${PAUSE_SECS:-600}
PROBE_URL=$(tail -n +2 "$RUNLIST" | head -1 | awk -F'\t' '{print $1}')
export BREAK_AFTER PAUSE_SECS PROBE_URL
: > logs/.consec_fail

note_fail() { ( flock 9; echo x >> logs/.consec_fail ) 9> logs/.cb.lock; }
note_ok()   { ( flock 9; : > logs/.consec_fail ) 9> logs/.cb.lock; }
wait_if_broken() {
  while [ "$(wc -l < logs/.consec_fail)" -ge "$BREAK_AFTER" ]; do
    (
      flock 9
      # re-check under the lock: another worker may already have cleared it
      [ "$(wc -l < logs/.consec_fail)" -ge "$BREAK_AFTER" ] || exit 0
      echo "[$(date '+%F %T')] circuit breaker OPEN: $(wc -l < logs/.consec_fail) consecutive failures; pausing ${PAUSE_SECS}s" >> logs/breaker.log
      sleep "$PAUSE_SECS"
      if curl -fs --max-time 180 -o /dev/null "$PROBE_URL" 2>/dev/null; then
        echo "[$(date '+%F %T')] probe OK, resuming" >> logs/breaker.log
        : > logs/.consec_fail
      else
        echo "[$(date '+%F %T')] probe still failing" >> logs/breaker.log
      fi
    ) 9> logs/.cb.lock
  done
}
export -f note_fail note_ok wait_if_broken

download_one() {
  sample=$1; url=$2; expected=$3
  [ -s quant/$sample/quant.sf ] && return 0
  [ -f fastq/$sample/.ready ] && return 0
  while [ "$(ls fastq 2>/dev/null | while read d; do [ -f fastq/$d/.ready ] && echo x; done | wc -l)" -ge "$BUFFER" ]; do sleep 20; done
  rm -rf fastq/$sample; mkdir -p fastq/$sample   # start clean: leftovers from an interrupted download must not survive
  tarfile="fastq/$sample/$sample.fastq.tar"
  ok=0
  for attempt in 1 2 3 4 5 6; do
    wait_if_broken
    curl -fsS --retry 2 --retry-delay 20 -o "$tarfile" "$url"
    rc=$?
    actual=$(stat -c%s "$tarfile" 2>/dev/null || echo 0)
    if [ "$rc" -eq 0 ] && { [ -z "$expected" ] || [ "$actual" = "$expected" ]; }; then ok=1; note_ok; break; fi
    note_fail
    echo "$sample attempt $attempt rc=$rc actual=$actual expected=$expected" >> logs/retries.log
    sleep $((attempt * 15))
  done
  if [ "$ok" != "1" ]; then
    rm -rf fastq/$sample; echo "$sample download-failed-after-retries" >> logs/failures.log; return 1
  fi
  # NeMO tars come in two layouts: some wrap the pair in a directory, some put the two
  # .fastq.gz files at the top level. --strip-components=1 silently extracted NOTHING from
  # the latter (it stripped the filename itself), so extract as-is and then find the pair.
  tar -xf "$tarfile" -C fastq/$sample || { rm -rf fastq/$sample; echo "$sample tar-extract-failed" >> logs/failures.log; return 1; }
  rm -f "$tarfile"
  r1=$(find fastq/$sample -type f -name '*_R1*.f*q.gz' | head -1)
  r2=$(find fastq/$sample -type f -name '*_R2*.f*q.gz' | head -1)
  [ -n "$r1" ] && [ "$r1" != "fastq/$sample/${sample}_R1.fastq.gz" ] && mv "$r1" fastq/$sample/${sample}_R1.fastq.gz
  [ -n "$r2" ] && [ "$r2" != "fastq/$sample/${sample}_R2.fastq.gz" ] && mv "$r2" fastq/$sample/${sample}_R2.fastq.gz
  if [ -s fastq/$sample/${sample}_R1.fastq.gz ] && [ -s fastq/$sample/${sample}_R2.fastq.gz ]; then
    touch fastq/$sample/.ready
  else
    rm -rf fastq/$sample; echo "$sample missing-R1R2" >> logs/failures.log
  fi
}
export -f download_one; export BUFFER
( tail -n +2 "$RUNLIST" | awk -F'\t' '{print $4"\t"$1"\t"$2}' \
    | xargs -P "$N_DL" -L 1 bash -c 'download_one "$0" "$1" "$2"'; touch fastq/.all_downloaded ) &
DL_PID=$!
log "downloader started (pid $DL_PID, $N_DL parallel)"

# ---------------------------------------------------------------- 3. quantifier loop
quant_one() {
  sample=$1
  r1=fastq/$sample/${sample}_R1.fastq.gz
  r2=fastq/$sample/${sample}_R2.fastq.gz
  salmon quant -i "$IDX" -l A -1 "$r1" -2 "$r2" -p "$T_Q" --validateMappings --seqBias --gcBias -o quant/$sample -q > logs/$sample.log 2>&1
  if [ -s quant/$sample/quant.sf ]; then
    rm -rf fastq/$sample; echo "$sample done $(grep -o '"percent_mapped": [0-9.]*' quant/$sample/aux_info/meta_info.json)" >> logs/progress.log
  else
    rm -rf fastq/$sample; rm -rf quant/$sample; echo "$sample quant-failed" >> logs/failures.log
  fi
}
export -f quant_one; export IDX T_Q
# Streaming scheduler: keep N_Q salmon jobs running at all times, starting a new one the moment
# any finishes. (The earlier batch-per-pass design waited for the slowest cell in each batch,
# which left most workers idle whenever cell depth is uneven — e.g. Tasic 2016.)
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
log "quantification finished: $(ls quant | wc -l) cells; failures: $(wc -l < logs/failures.log 2>/dev/null || echo 0)"

# ---------------------------------------------------------------- 4. collapse to tables
python3 - "$RUNLIST" <<'EOF'
import sys
import pandas as pd, glob, gzip, os, re, json
genes={'Kcnd1','Kcnd2','Kcnd3','Dpp6','Dpp10','Kcnip1','Kcnip2','Kcnip3','Kcnip4'}
tx={}
with gzip.open('ref/gencode.vM36.transcripts.fa.gz','rt') as f:
    for l in f:
        if l[0]=='>':
            p=l[1:].strip().split('|')
            if p[5] in genes: tx[p[0]]=(p[5],p[4],p[7] if len(p)>7 else '')
rows=[]; qc=[]
for q in sorted(glob.glob('quant/*/quant.sf')):
    sample=q.split('/')[1]; d=pd.read_csv(q,sep='\t',index_col=0)
    try:
        m=json.load(open(f'quant/{sample}/aux_info/meta_info.json')); qc.append((sample,m.get('num_processed'),m.get('percent_mapped')))
    except Exception: pass
    for t,(g,name,btype) in tx.items():
        if t in d.index: rows.append((sample,g,name,t,btype,d.loc[t,'TPM'],d.loc[t,'NumReads'],d.loc[t,'EffectiveLength']))
T=pd.DataFrame(rows,columns=['sample_name','gene','transcript','ensmust','biotype','TPM','reads','eff_len'])
pd.DataFrame(qc,columns=['sample_name','num_processed','percent_mapped']).to_csv('tables/run_qc.csv',index=False)
runlist=pd.read_csv(sys.argv[1],sep='\t')
J=T.merge(runlist[['sample_name','subclass_label','cluster_label','region_label']],on='sample_name',how='left')
print('matched cells:',J.sample_name.nunique(),'of',T.sample_name.nunique())
J['gene_tpm']=J.groupby(['sample_name','gene']).TPM.transform('sum')
J['frac']=J.TPM/J.gene_tpm.where(J.gene_tpm>0)
J.to_csv('tables/hip_isoform_per_cell_annotated.csv',index=False)
out=J.groupby(['gene','transcript','ensmust','subclass_label']).agg(
    mean_TPM=('TPM','mean'),mean_frac=('frac','mean'),pct_cells_detected=('reads',lambda x:(x>=3).mean()*100),n=('sample_name','nunique')).reset_index()
out.to_csv('tables/hip_isoform_by_subclass.csv',index=False)
print('wrote tables/')
EOF
log "all done"
