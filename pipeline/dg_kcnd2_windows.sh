#!/usr/bin/env bash
# Fine-scale coverage across the Kcnd2 locus in the DG single cells aligned by run_dg_star.sh.
# The Kv4-null cells carry ~94% of their intronic Kcnd2 depth in intron 1 (~400 kb). Two readings differ in WHERE in
# intron 1 the reads are: attenuated Kcnd2 transcription piles up just downstream of exon 1; an independent
# transcription unit inside the intron (lncRNA / eRNA / alternative promoter) makes a hotspot away from exon 1, present
# in null and high cells alike. This script bins the gene span (+/- 20 kb) into 10-kb windows and tabulates depth per
# cell, and lists every annotated gene/transcript overlapping the locus on either strand.
# Usage (in ~/kv4_dendrite):  bash dg_kcnd2_windows.sh   -> tables/dg_kcnd2_windows.tsv, tables/kcnd2_locus_annotation.tsv
set -uo pipefail
W=${W:-10000}
read -r CHR GS GE STRAND < <(awk -F'\t' '$4=="Kcnd2"{print $1,$2,$3,$6}' ref/kv4_genes.bed)
[ -n "${CHR:-}" ] || { echo "Kcnd2 not in ref/kv4_genes.bed"; exit 1; }
S=$((GS - 20000)); [ $S -lt 0 ] && S=0; E=$((GE + 20000))
echo "Kcnd2 locus $CHR:$S-$E ($STRAND), ${W}-bp windows"
awk -v c=$CHR -v s=$S -v e=$E -v w=$W 'BEGIN{OFS="\t"; i=0; for(x=s;x<e;x+=w){print c,x,(x+w<e?x+w:e),"w"i; i++}}' > ref/kcnd2_windows.bed
# exon / intron landmarks in the same coordinates, for the plot
{ awk -F'\t' 'BEGIN{OFS="\t"} $1=="Kcnd2"{print "exon",$3-1,$4}' <(tail -n +2 ref/kv4_exons.saf);
  awk -F'\t' 'BEGIN{OFS="\t"} $4 ~ /^Kcnd2_intron_/{print $4,$2,$3}' ref/kv4_introns.bed; } > tables/kcnd2_landmarks.tsv
# annotation overlapping the locus (both strands)
zcat ref/gencode.vM36.primary_assembly.annotation.gtf.gz | awk -F'\t' -v c=$CHR -v s=$S -v e=$E 'BEGIN{OFS="\t"} ($3=="gene"||$3=="transcript") && $1==c && $5>s && $4<e {
  match($9,/gene_name "[^"]+"/); gn=substr($9,RSTART+11,RLENGTH-12);
  match($9,/gene_type "[^"]+"/); gt=substr($9,RSTART+11,RLENGTH-12);
  tn=""; if(match($9,/transcript_name "[^"]+"/)) tn=substr($9,RSTART+17,RLENGTH-18);
  print $3,gn,gt,tn,$4,$5,$7}' > tables/kcnd2_locus_annotation.tsv
echo "annotation overlapping the locus:"; awk -F'\t' '$1=="gene"' tables/kcnd2_locus_annotation.tsv
# per-cell windowed depth
{ printf 'sample_name\twindow\tchr\tstart\tend\tdepth_sum\n'
  for d in star_dg/*/; do s=${d%/}; s=${s#star_dg/}
    [ -s star_dg/$s/region.bam ] || continue
    samtools bedcov ref/kcnd2_windows.bed star_dg/$s/region.bam | awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$4,$1,$2,$3,$5}'
  done; } > tables/dg_kcnd2_windows.tsv
echo "$(cut -f1 tables/dg_kcnd2_windows.tsv | sort -u | wc -l) cells -> tables/dg_kcnd2_windows.tsv"
