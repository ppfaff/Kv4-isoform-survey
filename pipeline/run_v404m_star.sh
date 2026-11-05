#!/usr/bin/env bash
# Genome alignment of the V404M vs wild-type bulk libraries: gene counts, splice junctions, Kv4 exon/intron
# coverage and a pileup at the V404M codon.
#
# Per sample it keeps only small outputs — ReadsPerGene.out.tab (an alignment-based count matrix, independent of
# salmon, for the differential expression), SJ.out.tab (genome-wide junctions: Xbp1 splicing for the ISR readout,
# the dFosB junction, Kv4.3L and the Kcnip3 acceptor), a BAM restricted to the Kv4-complex/Xbp1/Fosb loci, and
# per-exon/per-intron depth — then deletes the full alignment. The genome-wide BAM is never sorted (reads are
# filtered to the loci first), which is what keeps the run to roughly 20 minutes per sample.
#
# Prerequisites:  ref/star_idx (built by run_dendrite_star.sh), ref/kv4_*.saf|bed (make_intron_saf.py),
#                 tables/v404m_samples.tsv, tables/v404m_codon.tsv, and STAR/samtools/featureCounts on PATH.
# Usage:  T=16 nohup nice -n 19 bash run_v404m_star.sh > run_v404m_star.log 2>&1 &
#   N_A=2 runs two alignments at once (the genome is held in shared memory, so this costs no extra RAM); with 32
#   cores, N_A=2 T=12 is a reasonable setting. Everything is resumable — finished samples are skipped.
set -uo pipefail
T=${T:-16}
N_A=${N_A:-1}
BUFFER=${BUFFER:-2}
SAMPLES=${SAMPLES:-tables/v404m_samples.tsv}
STAGE_DIR=${STAGE_DIR:-fastq_star}
EXTRA_GENES=${EXTRA_GENES:-Xbp1 Fosb Atf4 Ddit3 Hspa5 Trib3 Chac1 Gfap}
ROOT=$(pwd)
mkdir -p "$STAGE_DIR" star_v404m logs tables; rm -f "$STAGE_DIR/.all_staged"
log(){ echo "[$(date '+%F %T')] $*"; }
for tool in STAR samtools featureCounts; do command -v $tool >/dev/null || { log "$tool not on PATH"; exit 1; }; done
[ -s ref/star_idx/SA ] || { log "ref/star_idx missing — build it with run_dendrite_star.sh"; exit 1; }
[ -s "$SAMPLES" ] || { log "missing $SAMPLES"; exit 1; }
GTF=${GTF:-ref/gencode.vM36.primary_assembly.annotation.gtf.gz}
[ -s ref/kv4_introns.saf ] || python3 "$(dirname "$0")/make_intron_saf.py" "$GTF"
[ -s ref/kv4_exons.bed ] || tail -n +2 ref/kv4_exons.saf | awk 'BEGIN{OFS="\t"}{print $2,$3-1,$4,$1"_exon_"NR,0,$5}' > ref/kv4_exons.bed

# regions kept in the per-sample BAM: the 9 Kv4-complex genes plus the stress genes we want coverage for
if [ ! -s ref/v404m_regions.bed ]; then
  cp ref/kv4_genes.bed ref/v404m_regions.bed
  for g in $EXTRA_GENES; do
    zcat "$GTF" | awk -F'\t' -v g="$g" 'BEGIN{OFS="\t"} $3=="gene" && $9 ~ "gene_name \""g"\"" {print $1,($4>2001?$4-2001:0),$5+2000,g,0,$7}' >> ref/v404m_regions.bed
  done
  sort -k1,1 -k2,2n ref/v404m_regions.bed -o ref/v404m_regions.bed
  log "regions: $(wc -l < ref/v404m_regions.bed) loci in ref/v404m_regions.bed"
fi
TOTAL=$(tail -n +2 "$SAMPLES" | wc -l); log "$TOTAL samples"

# ---------------------------------------------------------------- stager (background)
stage_one() {
  s=$1; r1=$2; r2=$3
  [ -s "star_v404m/$s/ReadsPerGene.out.tab" ] && return 0
  [ -f "$STAGE_DIR/$s/.ready" ] && return 0
  while [ "$(find "$STAGE_DIR" -name .ready 2>/dev/null | wc -l)" -ge "$BUFFER" ]; do sleep 20; done
  rm -rf "${STAGE_DIR:?}/$s"; mkdir -p "$STAGE_DIR/$s"
  if cp "$r1" "$STAGE_DIR/$s/R1.fastq.gz" && cp "$r2" "$STAGE_DIR/$s/R2.fastq.gz"; then touch "$STAGE_DIR/$s/.ready"
  else rm -rf "${STAGE_DIR:?}/$s"; echo "$s stage-failed" >> logs/v404m_star_failures.log; fi
}
(
  tail -n +2 "$SAMPLES" | while IFS=$'\t' read -r s tissue genotype sex r1 r2 rest; do
    [ -n "$s" ] || continue; stage_one "$s" "$r1" "$r2"
  done
  touch "$STAGE_DIR/.all_staged"
) &
log "stager started (pid $!)"

# ---------------------------------------------------------------- genome into shared memory
STAR --genomeLoad LoadAndExit --genomeDir ref/star_idx --outFileNamePrefix logs/_gload. > logs/v404m_genome_load.log 2>&1 \
  || { log "genome load failed (see logs/v404m_genome_load.log)"; exit 1; }
log "genome loaded into shared memory"
cleanup(){ STAR --genomeLoad Remove --genomeDir ref/star_idx --outFileNamePrefix logs/_gunload. >/dev/null 2>&1; }
trap cleanup EXIT

CHROM=$(awk -F'\t' '$1=="chrom"{print $2}' tables/v404m_codon.tsv 2>/dev/null)

align_one() {
  s=$1
  d=star_v404m/$s; mkdir -p "$d"
  STAR --runThreadN "$T" --genomeDir ref/star_idx --genomeLoad LoadAndKeep \
       --readFilesIn "$STAGE_DIR/$s/R1.fastq.gz" "$STAGE_DIR/$s/R2.fastq.gz" --readFilesCommand zcat \
       --outFileNamePrefix "$d/" --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --quantMode GeneCounts \
       --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
       --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
       > "logs/$s.star.log" 2>&1
  bam=$d/Aligned.out.bam
  if [ ! -s "$bam" ]; then echo "$s star-failed" >> logs/v404m_star_failures.log; rm -rf "${STAGE_DIR:?}/$s"; return; fi
  # restrict to the loci of interest BEFORE sorting — the genome-wide sort is the expensive step and is not needed
  samtools view -@ 2 -b -L ref/v404m_regions.bed -o "$d/region.unsorted.bam" "$bam" \
    && samtools sort -@ 2 -m 2G -o "$d/region.bam" "$d/region.unsorted.bam" && samtools index "$d/region.bam"
  rm -f "$bam" "$d/region.unsorted.bam"
  for st in 0 1 2; do
    featureCounts -T 2 -p --countReadPairs -s $st -F SAF -a ref/kv4_exons.saf   -o "$d/exons.s$st.txt"   "$d/region.bam" > "logs/$s.fc_e$st.log" 2>&1
    featureCounts -T 2 -p --countReadPairs -s $st -F SAF -a ref/kv4_introns.saf -o "$d/introns.s$st.txt" "$d/region.bam" > "logs/$s.fc_i$st.log" 2>&1
  done
  samtools bedcov ref/kv4_introns.bed "$d/region.bam" > "$d/intron_cov.tsv"
  samtools bedcov ref/kv4_exons.bed   "$d/region.bam" > "$d/exon_cov.tsv"
  if [ -n "$CHROM" ]; then
    samtools mpileup -f ref/GRCm39.primary_assembly.genome.fa -l tables/v404m_codon.pos -d 100000 -Q 13 \
      "$d/region.bam" > "$d/v404m_pileup.txt" 2>/dev/null
  fi
  rm -rf "$d/_STARtmp" "${STAGE_DIR:?}/$s"
  echo "$s done $(grep 'Uniquely mapped reads %' "$d/Log.final.out" | tr -s ' ' | cut -d'|' -f2)" >> logs/v404m_star_progress.log
}
export -f align_one; export T STAGE_DIR CHROM

last=0
while true; do
  running=$(jobs -rp | wc -l)
  if [ "$running" -lt "$N_A" ]; then
    next=$(for d in "$STAGE_DIR"/*/; do d=${d%/}; [ -f "$d/.ready" ] && { basename "$d"; break; }; done)
    if [ -n "$next" ]; then mv "$STAGE_DIR/$next/.ready" "$STAGE_DIR/$next/.claimed"; align_one "$next" & continue; fi
  fi
  if [ "$running" -eq 0 ] && [ -f "$STAGE_DIR/.all_staged" ] && [ -z "$(find "$STAGE_DIR" -name .ready -print -quit 2>/dev/null)" ]; then break; fi
  now=$(date +%s)
  if [ $((now - last)) -ge 600 ]; then
    log "aligned $(ls star_v404m 2>/dev/null | wc -l) / $TOTAL (running $running, staged $(find "$STAGE_DIR" -name .ready 2>/dev/null | wc -l))"
    last=$now
  fi
  sleep 20
done
wait

# ---------------------------------------------------------------- collect the small per-sample outputs
log "collecting"
python3 - "$SAMPLES" <<'EOF'
import sys, glob, os, re
import pandas as pd
S = pd.read_csv(sys.argv[1], sep='\t')
meta = S.set_index('sample')[['tissue', 'genotype', 'sex']]
# gene counts (STAR column 2 = unstranded, 3 = fwd, 4 = rev)
counts, strand_qc = {}, []
for f in glob.glob('star_v404m/*/ReadsPerGene.out.tab'):
    s = f.split('/')[1]
    t = pd.read_csv(f, sep='\t', header=None, names=['gene_id', 'unstr', 'fwd', 'rev'])
    t = t[~t.gene_id.str.startswith('N_')]
    tot = t[['unstr', 'fwd', 'rev']].sum()
    strand_qc.append((s, *tot.values, tot['fwd'] / max(tot['fwd'] + tot['rev'], 1)))
    counts[s] = t.set_index('gene_id')
QC = pd.DataFrame(strand_qc, columns=['sample', 'assigned_unstranded', 'assigned_fwd', 'assigned_rev', 'frac_fwd']).set_index('sample').join(meta)
QC.to_csv('tables/v404m_star_qc.csv')
print(QC.round(3).to_string())
col = 'rev' if QC.frac_fwd.median() < 0.3 else ('fwd' if QC.frac_fwd.median() > 0.7 else 'unstr')
print(f'strand column used for the count matrix: {col} (median fwd fraction {QC.frac_fwd.median():.3f})')
if counts:
    M = pd.DataFrame({s: d[col] for s, d in counts.items()})
    M.index.name = 'gene_id'
    M.to_csv('tables/v404m_star_counts.csv.gz', compression='gzip')
    print(f'count matrix {M.shape[0]} genes x {M.shape[1]} samples -> tables/v404m_star_counts.csv.gz')
# junctions, exon/intron, pileup
rows = []
for f in glob.glob('star_v404m/*/SJ.out.tab'):
    s = f.split('/')[1]
    if os.path.getsize(f) == 0: continue
    t = pd.read_csv(f, sep='\t', header=None, names=['chr', 'start', 'end', 'strand', 'motif', 'annot', 'uniq', 'multi', 'overhang'])
    t['sample'] = s; rows.append(t)
if rows:
    J = pd.concat(rows)
    J.to_csv('tables/v404m_junctions_all.tsv.gz', sep='\t', index=False, compression='gzip')
    print(f'junctions: {len(J):,} rows -> tables/v404m_junctions_all.tsv.gz')
ei = []
for f in glob.glob('star_v404m/*/exons.s*.txt') + glob.glob('star_v404m/*/introns.s*.txt'):
    s = f.split('/')[1]; kind, st = re.match(r'(exons|introns)\.s(\d)\.txt', os.path.basename(f)).groups()
    t = pd.read_csv(f, sep='\t', comment='#'); t = t.rename(columns={t.columns[-1]: 'count'})
    g = t.groupby(t.Geneid.str.replace(r'_intron_\d+$', '', regex=True))['count'].sum()
    for gn, c in g.items(): ei.append((s, gn, kind[:-1], int(st), int(c)))
if ei:
    pd.DataFrame(ei, columns=['sample', 'gene', 'feature', 'strand_setting', 'count']).to_csv('tables/v404m_exon_intron.tsv', sep='\t', index=False)
    print('exon/intron -> tables/v404m_exon_intron.tsv')
pu = []
for f in glob.glob('star_v404m/*/v404m_pileup.txt'):
    s = f.split('/')[1]
    for line in open(f):
        p = line.rstrip('\n').split('\t')
        if len(p) >= 5: pu.append((s, p[0], int(p[1]), p[2], int(p[3]), p[4]))
if pu:
    pd.DataFrame(pu, columns=['sample', 'chr', 'pos', 'ref', 'depth', 'bases']).to_csv('tables/v404m_pileup.tsv', sep='\t', index=False)
    print('pileup -> tables/v404m_pileup.tsv')
EOF
cp ref/v404m_regions.bed ref/kv4_introns.bed ref/kv4_exons.bed tables/ 2>/dev/null
log "done — tables/v404m_star_counts.csv.gz, v404m_junctions_all.tsv.gz, v404m_exon_intron.tsv, v404m_pileup.tsv, v404m_star_qc.csv"
