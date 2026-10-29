#!/usr/bin/env bash
# Genome alignment of the Farris 2019 laminar libraries (cell-body vs dendritic layers) for
# exonic-vs-intronic Kv4-complex reads and per-intron Kcnd2 coverage by compartment.
#
# Steps (all resumable):
#   1. reference: GRCm39 primary assembly + GENCODE vM36 GTF -> STAR index (needs ~32 GB RAM, ~1 h on 16 threads)
#   2. per run: FASTQ (re-used from fastq/<run> if run_dendrite_isoforms.sh was run with KEEP_FASTQ=1,
#      otherwise downloaded with prefetch/fasterq-dump) -> STAR sorted BAM -> featureCounts on merged-exon
#      and intron SAF (all three strand settings) -> gene-region BAM kept, full BAM deleted
#   3. per-intron depth (samtools bedcov) on the region BAMs -> tables/star_intron_depth.tsv
#   4. summary tables: tables/star_exon_intron_counts.tsv, tables/star_ReadsPerGene.tsv
#
# Dependencies (into the kv4 env):  micromamba install -c bioconda -c conda-forge star=2.7.11b subread samtools
# Usage:  nohup nice -n 19 bash run_dendrite_star.sh > run_star.log 2>&1 &
#   RUNS=tables/farris_runs.txt (default: all 'farris' rows of tables/dendrite_samples.tsv)
#   T=16 threads for STAR/featureCounts; STAR_RAM=48000000000 for index generation
set -uo pipefail
T=${T:-16}
STAR_RAM=${STAR_RAM:-48000000000}
SAMPLES=${SAMPLES:-tables/dendrite_samples.tsv}
RUNS=${RUNS:-}
ROOT=$(pwd)
mkdir -p ref star fastq sra logs tables
log(){ echo "[$(date '+%F %T')] $*"; }
for tool in STAR featureCounts samtools; do command -v $tool >/dev/null || { log "$tool not on PATH"; exit 1; }; done

if [ -z "$RUNS" ]; then
  RUNS=tables/farris_runs.txt
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="dataset")d=i;if($i=="run_accession")r=i}} NR>1&&$d=="farris"{print $r}' "$SAMPLES" > "$RUNS"
fi
[ -s "$RUNS" ] || { log "no runs in $RUNS"; exit 1; }
TOTAL=$(wc -l < "$RUNS"); log "$TOTAL Farris runs"

# ---------------------------------------------------------------- 1. reference + STAR index
GC=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36
cd ref
[ -s GRCm39.primary_assembly.genome.fa ] || { curl -sS -O $GC/GRCm39.primary_assembly.genome.fa.gz && gunzip -f GRCm39.primary_assembly.genome.fa.gz; }
[ -s gencode.vM36.primary_assembly.annotation.gtf.gz ] || curl -sS -O $GC/gencode.vM36.primary_assembly.annotation.gtf.gz
[ -s gencode.vM36.primary_assembly.annotation.gtf ] || gunzip -k gencode.vM36.primary_assembly.annotation.gtf.gz
if [ ! -s star_idx/SA ]; then
  log "building STAR index (this is the slow step)"
  mkdir -p star_idx
  STAR --runMode genomeGenerate --runThreadN "$T" --genomeDir star_idx \
       --genomeFastaFiles GRCm39.primary_assembly.genome.fa \
       --sjdbGTFfile gencode.vM36.primary_assembly.annotation.gtf --sjdbOverhang 100 \
       --limitGenomeGenerateRAM "$STAR_RAM" > ../logs/star_index.log 2>&1 || { log "STAR index failed, see logs/star_index.log"; exit 1; }
fi
cd "$ROOT"
[ -s ref/kv4_introns.saf ] || python3 "$(dirname "$0")/make_intron_saf.py" ref/gencode.vM36.primary_assembly.annotation.gtf.gz

# ---------------------------------------------------------------- 2. download helper (same as the salmon script)
get_fastq() {   # $1 = run; leaves fastq/$run/<run>_1.fastq, _2.fastq
  run=$1
  if ls fastq/$run/${run}_1.f*q* >/dev/null 2>&1 && ls fastq/$run/${run}_2.f*q* >/dev/null 2>&1; then return 0; fi
  rm -rf fastq/$run; mkdir -p fastq/$run
  for attempt in 1 2 3; do
    timeout -k 30 5400 prefetch -q --max-size u -O sra "$run" >> logs/$run.prefetch.log 2>&1 && break
    rm -rf sra/$run; sleep 60
  done
  fasterq-dump -q -e 4 --split-3 -t sra/$run/tmp -O fastq/$run sra/$run > logs/$run.dump.log 2>&1
  rm -rf sra/$run
  ls fastq/$run/${run}_1.fastq fastq/$run/${run}_2.fastq >/dev/null 2>&1
}

# ---------------------------------------------------------------- 3. align + count, one run at a time (STAR uses all threads)
while read -r run; do
  [ -n "$run" ] || continue
  if [ -s star/$run/region.bam.bai ] && [ -s star/$run/exons.s2.txt ]; then continue; fi
  log "$run: fetching FASTQ"
  get_fastq "$run" || { echo "$run fastq-failed" >> logs/star_failures.log; log "$run: FASTQ failed"; continue; }
  r1=$(ls fastq/$run/${run}_1.f*q*); r2=$(ls fastq/$run/${run}_2.f*q*)
  case "$r1" in *.gz) RC="--readFilesCommand zcat";; *) RC="";; esac
  mkdir -p star/$run
  log "$run: STAR"
  STAR --runThreadN "$T" --genomeDir ref/star_idx --readFilesIn "$r1" "$r2" $RC \
       --outFileNamePrefix star/$run/ --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif \
       --quantMode GeneCounts --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
       --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
       --limitBAMsortRAM 20000000000 > logs/$run.star.log 2>&1
  bam=star/$run/Aligned.sortedByCoord.out.bam
  [ -s "$bam" ] || { echo "$run star-failed" >> logs/star_failures.log; log "$run: STAR failed"; continue; }
  samtools index -@ 4 "$bam"
  # exon / intron counts for the 9 genes, all strand settings so strandedness can be verified from the counts
  for s in 0 1 2; do
    featureCounts -T 4 -p --countReadPairs -s $s -F SAF -a ref/kv4_exons.saf   -o star/$run/exons.s$s.txt   "$bam" > logs/$run.fc_exons.s$s.log 2>&1
    featureCounts -T 4 -p --countReadPairs -s $s -F SAF -a ref/kv4_introns.saf -o star/$run/introns.s$s.txt "$bam" > logs/$run.fc_introns.s$s.log 2>&1
  done
  # keep only the reads in the 9 gene regions (+/- 2 kb), drop the full BAM and the FASTQ
  samtools view -@ 4 -b -L ref/kv4_genes.bed -o star/$run/region.bam "$bam" && samtools index star/$run/region.bam
  rm -f "$bam" "$bam.bai"; rm -rf fastq/$run star/$run/_STARtmp
  log "$run: done  $(grep 'Uniquely mapped reads %' star/$run/Log.final.out | tr -s ' ' | cut -d'|' -f2)"
done < "$RUNS"

# ---------------------------------------------------------------- 4. summaries
log "per-intron depth"
{ printf 'run\tintron\tchr\tstart\tend\tdepth_sum\tmean_depth\n'
  for d in star/*/; do run=${d%/}; run=${run#star/}
    [ -s star/$run/region.bam ] || continue
    samtools bedcov ref/kv4_introns.bed star/$run/region.bam | awk -v r="$run" 'BEGIN{OFS="\t"}{print r,$4,$1,$2,$3,$7,$7/($3-$2)}'
  done; } > tables/star_intron_depth.tsv

python3 - <<'EOF'
import glob, pandas as pd, os, re
rows = []
for f in glob.glob('star/*/exons.s*.txt') + glob.glob('star/*/introns.s*.txt'):
    run = f.split('/')[1]; kind, s = re.match(r'(exons|introns)\.s(\d)\.txt', os.path.basename(f)).groups()
    t = pd.read_csv(f, sep='\t', comment='#')
    t = t.rename(columns={t.columns[-1]: 'count'})
    gene = t.Geneid.str.replace(r'_intron_\d+$', '', regex=True)
    g = t.groupby(gene)['count'].sum()
    for gn, c in g.items():
        rows.append((run, gn, kind[:-1], int(s), int(c)))
E = pd.DataFrame(rows, columns=['run', 'gene', 'feature', 'strand_setting', 'count'])
E.to_csv('tables/star_exon_intron_counts.tsv', sep='\t', index=False)
# STAR's ReadsPerGene (unstranded, stranded-forward, stranded-reverse) for every gene -> strandedness check + library size
R = []
for f in glob.glob('star/*/ReadsPerGene.out.tab'):
    run = f.split('/')[1]
    t = pd.read_csv(f, sep='\t', header=None, names=['gene_id', 'unstranded', 'fwd', 'rev'])
    tot = t[~t.gene_id.str.startswith('N_')][['unstranded', 'fwd', 'rev']].sum()
    R.append((run, *tot.values))
pd.DataFrame(R, columns=['run', 'assigned_unstranded', 'assigned_fwd', 'assigned_rev']).to_csv('tables/star_ReadsPerGene.tsv', sep='\t', index=False)
print(E.pivot_table(index=['run', 'gene'], columns=['feature', 'strand_setting'], values='count').head(20))
EOF
log "all done — tables/star_exon_intron_counts.tsv, tables/star_intron_depth.tsv, tables/star_ReadsPerGene.tsv; run scripts/dendrite_compartment.py"
