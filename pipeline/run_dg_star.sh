#!/usr/bin/env bash
# Read-level Kcnd2 analysis in single DG granule cells (Allen SSv4, Yao 2021): STAR-align ~750 cells (all Kv4-null
# cells + Kv4-low and Kv4-high controls), keep only reads in the 9 Kv4-complex gene regions, and tabulate per cell
#   * per-intron and per-exon-block depth for Kcnd2 / Kcnip4 / Dpp6  (samtools bedcov)
#   * splice-junction counts in those genes                           (SJ.out.tab)
# so that the Kv4-null cells can be classified: nascent transcription caught mid-gene (5'-biased intronic reads),
# a detained intron (one intron high, others spliced), or a post-splicing loss (flat introns, residual junction reads).
#
# Needs the STAR index from run_dendrite_star.sh (ref/star_idx) and the SAF/BED files from make_intron_saf.py; run in
# the same work dir (~/kv4_dendrite) so they are re-used. The index is loaded ONCE into shared memory (STAR
# --genomeLoad LoadAndKeep), which is what makes 750 small alignments feasible (~20-40 s per cell instead of 2 min).
#
# Usage:  RUNLIST=tables/dg_readlevel_runlist.tsv T=16 nohup nice -n 19 bash run_dg_star.sh > run_dg_star.log 2>&1 &
set -uo pipefail
T=${T:-16}
N_DL=${N_DL:-2}
RUNLIST=${RUNLIST:-tables/dg_readlevel_runlist.tsv}
ROOT=$(pwd)
mkdir -p fastq_dg star_dg logs tables
log(){ echo "[$(date '+%F %T')] $*"; }
for tool in STAR samtools; do command -v $tool >/dev/null || { log "$tool not on PATH"; exit 1; }; done
[ -s ref/star_idx/SA ] || { log "ref/star_idx missing — run run_dendrite_star.sh first (index step)"; exit 1; }
[ -s ref/kv4_genes.bed ] || python3 "$(dirname "$0")/make_intron_saf.py" ref/gencode.vM36.primary_assembly.annotation.gtf.gz
[ -s ref/kv4_exons.bed ] || tail -n +2 ref/kv4_exons.saf | awk 'BEGIN{OFS="\t"}{print $2,$3-1,$4,$1"_exon_"NR,0,$5}' > ref/kv4_exons.bed
[ -s "$RUNLIST" ] || { log "missing $RUNLIST — run build_dg_runlist.py"; exit 1; }
TOTAL=$(tail -n +2 "$RUNLIST" | wc -l); log "$TOTAL cells"
H=$(head -1 "$RUNLIST"); col(){ echo "$H" | tr '\t' '\n' | grep -nx "$1" | cut -d: -f1; }
C_S=$(col sample_name); C_U=$(col fastq_url 2>/dev/null || true)
# the hip runlist's URL column name varies; fall back to the first column that looks like a URL
if [ -z "$C_U" ]; then C_U=$(tail -n +2 "$RUNLIST" | head -1 | tr '\t' '\n' | grep -n '^http' | head -1 | cut -d: -f1); fi
log "sample column $C_S, url column $C_U"

# ---------------------------------------------------------------- downloader (NeMO tar -> R1/R2), bounded buffer
download_one() {   # $1 = sample, $2 = url
  s=$1; url=$2
  [ -s star_dg/$s/region.bam ] && return 0
  [ -f fastq_dg/$s/.ready ] && return 0
  while [ "$(find fastq_dg -name .ready 2>/dev/null | wc -l)" -ge 20 ]; do sleep 10; done
  rm -rf fastq_dg/$s; mkdir -p fastq_dg/$s
  for attempt in 1 2 3 4; do
    curl -fsS --retry 2 --retry-delay 15 -o fastq_dg/$s/$s.tar "$url" && break
    sleep $((attempt * 20))
  done
  tar -xf fastq_dg/$s/$s.tar -C fastq_dg/$s 2>/dev/null; rm -f fastq_dg/$s/$s.tar
  r1=$(find fastq_dg/$s -type f -name '*_R1*.f*q.gz' | head -1); r2=$(find fastq_dg/$s -type f -name '*_R2*.f*q.gz' | head -1)
  if [ -n "$r1" ] && [ -n "$r2" ]; then mv "$r1" fastq_dg/$s/R1.fastq.gz; mv "$r2" fastq_dg/$s/R2.fastq.gz; touch fastq_dg/$s/.ready
  else rm -rf fastq_dg/$s; echo "$s download-failed" >> logs/dg_failures.log; fi
}
export -f download_one
( tail -n +2 "$RUNLIST" | awk -F'\t' -v s=$C_S -v u=$C_U '{print $s"\t"$u}' | xargs -P "$N_DL" -L 1 bash -c 'download_one "$0" "$1"'; touch fastq_dg/.all_downloaded ) &
log "downloader started (pid $!)"

# ---------------------------------------------------------------- STAR with the genome held in shared memory
STAR --genomeLoad LoadAndExit --genomeDir ref/star_idx > logs/dg_genome_load.log 2>&1 || { log "genome load failed (see logs/dg_genome_load.log; is enough RAM free?)"; exit 1; }
log "genome loaded into shared memory"
cleanup(){ STAR --genomeLoad Remove --genomeDir ref/star_idx > /dev/null 2>&1; }
trap cleanup EXIT

align_one() {   # $1 = sample
  s=$1
  mkdir -p star_dg/$s
  STAR --runThreadN "$T" --genomeDir ref/star_idx --genomeLoad LoadAndKeep \
       --readFilesIn fastq_dg/$s/R1.fastq.gz fastq_dg/$s/R2.fastq.gz --readFilesCommand zcat \
       --outFileNamePrefix star_dg/$s/ --outSAMtype BAM Unsorted --outSAMstrandField intronMotif \
       --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNoverReadLmax 0.04 \
       --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 > logs/$s.dgstar.log 2>&1
  bam=star_dg/$s/Aligned.out.bam
  [ -s "$bam" ] || { echo "$s star-failed" >> logs/dg_failures.log; rm -rf fastq_dg/$s; return; }
  samtools sort -@ 4 -o star_dg/$s/sorted.bam "$bam" 2>/dev/null && samtools index star_dg/$s/sorted.bam
  samtools view -b -L ref/kv4_genes.bed -o star_dg/$s/region.bam star_dg/$s/sorted.bam && samtools index star_dg/$s/region.bam
  samtools bedcov ref/kv4_introns.bed star_dg/$s/region.bam > star_dg/$s/intron_cov.tsv
  samtools bedcov ref/kv4_exons.bed   star_dg/$s/region.bam > star_dg/$s/exon_cov.tsv
  # junctions inside any of the 9 gene regions (intervals from kv4_genes.bed)
  awk -F'\t' 'NR==FNR{c[NR]=$1;a[NR]=$2;b[NR]=$3;n=NR;next} {for(i=1;i<=n;i++) if($1==c[i] && $2>=a[i] && $3<=b[i]){print;break}}' ref/kv4_genes.bed star_dg/$s/SJ.out.tab > star_dg/$s/sj_region.tsv
  rm -f "$bam" star_dg/$s/sorted.bam star_dg/$s/sorted.bam.bai; rm -rf star_dg/$s/_STARtmp fastq_dg/$s
  echo "$s done $(grep 'Uniquely mapped reads %' star_dg/$s/Log.final.out | tr -s ' ' | cut -d'|' -f2)" >> logs/dg_progress.log
}

while true; do
  next=$(for d in fastq_dg/*/; do d=${d%/}; [ -f "$d/.ready" ] && { echo "${d#fastq_dg/}"; break; }; done)
  if [ -n "$next" ]; then mv fastq_dg/$next/.ready fastq_dg/$next/.claimed; align_one "$next"; continue; fi
  [ -f fastq_dg/.all_downloaded ] && [ -z "$(find fastq_dg -name .ready -print -quit 2>/dev/null)" ] && break
  sleep 5
done
log "aligned: $(ls star_dg | wc -l) / $TOTAL; failures: $(wc -l < logs/dg_failures.log 2>/dev/null || echo 0)"

# ---------------------------------------------------------------- collect
python3 - "$RUNLIST" <<'EOF'
import sys, glob, os, pandas as pd
R = pd.read_csv(sys.argv[1], sep='\t', dtype=str)
grp = dict(zip(R.sample_name, R.grp))
cov, sj = [], []
for d in glob.glob('star_dg/*/'):
    s = d.strip('/').split('/')[-1]
    for kind in ['intron', 'exon']:
        f = f'{d}{kind}_cov.tsv'
        if os.path.exists(f):
            t = pd.read_csv(f, sep='\t', header=None, names=['chr', 'start', 'end', 'name', 'score', 'strand', 'depth_sum'])
            t['sample_name'] = s; t['grp'] = grp.get(s); t['feature'] = kind
            cov.append(t)
    f = f'{d}sj_region.tsv'
    if os.path.exists(f) and os.path.getsize(f):
        t = pd.read_csv(f, sep='\t', header=None, names=['chr', 'start', 'end', 'strand', 'motif', 'annot', 'uniq', 'multi', 'overhang'])
        t['sample_name'] = s; t['grp'] = grp.get(s); sj.append(t)
pd.concat(cov).to_csv('tables/dg_readlevel_coverage.tsv', sep='\t', index=False)
pd.concat(sj).to_csv('tables/dg_readlevel_junctions.tsv', sep='\t', index=False) if sj else None
print(f'{len(set(c.sample_name.iloc[0] for c in cov))} cells collected -> tables/dg_readlevel_coverage.tsv, tables/dg_readlevel_junctions.tsv')
EOF
log "done — run scripts/dg_readlevel.py on tables/dg_readlevel_coverage.tsv + dg_readlevel_junctions.tsv"
