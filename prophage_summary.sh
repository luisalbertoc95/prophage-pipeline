#!/bin/bash

# Summary of prophage predictions by source (GeNomad vs PhiSpy) and binning status

for dir in NovaSeq_*/; do
  echo "=== $dir ==="

  # MAG prophages (binned)
  echo -n "  MAG - GeNomad raw: "
  wc -l < "$dir/phage_analysis/mags/genomad_all.bed" 2>/dev/null || echo "0"
  echo -n "  MAG - PhiSpy all: "
  wc -l < "$dir/phage_analysis/mags/phispy_all.bed" 2>/dev/null || echo "0"
  echo -n "  MAG - PhiSpy unique: "
  wc -l < "$dir/phage_analysis/mags/phispy_unique.bed" 2>/dev/null || echo "0"
  echo -n "  MAG - merged total: "
  tail -n +2 "$dir/phage_analysis/mags/prophage_table.tsv" 2>/dev/null | wc -l

  # Unbinned prophages
  echo -n "  Unbinned - total: "
  tail -n +2 "$dir/phage_analysis/unbinned/prophage_table.tsv" 2>/dev/null | wc -l

  echo -n "  FINAL TOTAL: "
  tail -n +2 "$dir/phage_analysis/final_prophage_table.tsv" 2>/dev/null | wc -l
done
