#!/usr/bin/env bash
# Allele-specific Kcnd2 expression at the V404M codon, straight from the FASTQs (no alignment needed).
#
# For each sample it streams R1+R2 through a fixed-string grep for the four anchor patterns and tabulates the three
# bases sitting at codon 404. The wild-type and knock-in alleles are then simply two different triplets, and their
# read counts give the mutant fraction of Kcnd2 mRNA per sample — which says whether the mice are heterozygous or
# homozygous and whether the mutant allele is under-represented (e.g. by nonsense-mediated-decay-like effects,
# which a missense change should not cause).
#
# Prerequisites:  python3 find_v404m_codon.py        # -> tables/v404m_codon.tsv
#                 tables/v404m_samples.tsv
# Usage:  P=4 nohup nice -n 19 bash count_v404m_allele.sh > run_v404m_allele.log 2>&1 &
#   P=4 runs four samples at once; each is limited by gunzip throughput (~3 min per sample), so the full set takes
#   roughly 40 minutes at P=4. Run it after (or at low priority alongside) the salmon pass to avoid competing for
#   the network drive.
set -uo pipefail
P=${P:-4}
SAMPLES=${SAMPLES:-tables/v404m_samples.tsv}
CODON=${CODON:-tables/v404m_codon.tsv}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p logs tables allele
log(){ echo "[$(date '+%F %T')] $*"; }
[ -s "$CODON" ] || { log "missing $CODON — run find_v404m_codon.py first"; exit 1; }
[ -s "$SAMPLES" ] || { log "missing $SAMPLES"; exit 1; }
if command -v pigz >/dev/null 2>&1; then DECOMP="pigz -dc -p 2"; else DECOMP="zcat"; fi
log "decompressor: $DECOMP"

L=$(awk -F'\t' '$1=="left_anchor"{print $2}' "$CODON")
R=$(awk -F'\t' '$1=="right_anchor"{print $2}' "$CODON")
LRC=$(awk -F'\t' '$1=="left_anchor_rc"{print $2}' "$CODON")
RRC=$(awk -F'\t' '$1=="right_anchor_rc"{print $2}' "$CODON")
[ -n "$L" ] && [ -n "$R" ] || { log "anchors missing from $CODON"; exit 1; }
log "anchors: L=$L R=$R"
export L R LRC RRC DECOMP HERE CODON

one() {   # $1 sample  $2 r1  $3 r2
  s=$1; r1=$2; r2=$3
  [ -s "allele/$s.tsv" ] && return 0
  [ -s "$r1" ] && [ -s "$r2" ] || { echo "$s missing-fastq" >> logs/v404m_allele_failures.log; return 0; }
  $DECOMP "$r1" "$r2" 2>/dev/null \
    | LC_ALL=C grep -F -e "$L" -e "$R" -e "$LRC" -e "$RRC" \
    | python3 -I "$HERE/v404m_allele_extract.py" "$s" "$CODON" > "allele/$s.tsv.part" 2>> "logs/$s.allele.log"
  mv "allele/$s.tsv.part" "allele/$s.tsv"
  echo "[$(date '+%F %T')] $s $(awk '{n+=$5} END{print n+0}' "allele/$s.tsv") anchored reads" >> logs/v404m_allele_progress.log
}
export -f one

tail -n +2 "$SAMPLES" | awk -F'\t' '{print $1"\t"$5"\t"$6}' \
  | xargs -P "$P" -d '\n' -I{} bash -c 'IFS=$'"'"'\t'"'"' read -r s r1 r2 <<< "{}"; one "$s" "$r1" "$r2"'

{ printf 'sample\tanchor\torientation\ttriplet\tcount\twindow\n'; cat allele/*.tsv; } > tables/v404m_allele_counts.tsv
log "done: $(ls allele/*.tsv 2>/dev/null | wc -l) samples -> tables/v404m_allele_counts.tsv"
log "triplets seen overall:"
awk -F'\t' 'NR>1{n[$4]+=$5} END{for(t in n) print "   ", t, n[t]}' tables/v404m_allele_counts.tsv | sort -k3,3nr
