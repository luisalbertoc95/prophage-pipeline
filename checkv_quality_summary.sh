#!/bin/bash

# CheckV Quality Summary for All Prophages and Free Phages

echo "==============================================="
echo "CheckV Quality Assessment Summary"
echo "==============================================="
echo ""

# Per-sample quality distribution
for dir in NovaSeq_*/; do
  echo "=== $(basename $dir) ==="

  # All prophages (MAG + unbinned combined)
  if [ -f "$dir/phage_analysis/checkv_all_prophages/quality_summary.tsv" ]; then
    echo "ALL PROPHAGES (MAG + Unbinned):"
    echo "  Quality distribution:"
    tail -n +2 "$dir/phage_analysis/checkv_all_prophages/quality_summary.tsv" | cut -f8 | sort | uniq -c | \
      awk '{printf "    %-20s: %s\n", $2, $1}'

    echo -n "  Total assessed: "
    tail -n +2 "$dir/phage_analysis/checkv_all_prophages/quality_summary.tsv" | wc -l

    echo -n "  Mean completeness: "
    tail -n +2 "$dir/phage_analysis/checkv_all_prophages/quality_summary.tsv" | \
      awk -F'\t' '{if ($9 != "NA") sum+=$9; count++} END {if (count>0) printf "%.1f%%\n", sum/count; else print "N/A"}'

    echo -n "  Mean contamination: "
    tail -n +2 "$dir/phage_analysis/checkv_all_prophages/quality_summary.tsv" | \
      awk -F'\t' '{if ($10 != "NA") sum+=$10; count++} END {if (count>0) printf "%.1f%%\n", sum/count; else print "N/A"}'
  else
    echo "ALL PROPHAGES: CheckV results not found"
  fi

  # Free phages
  if [ -f "$dir/phage_analysis/checkv_free_phages/quality_summary.tsv" ]; then
    echo ""
    echo "FREE PHAGES:"
    echo "  Quality distribution:"
    tail -n +2 "$dir/phage_analysis/checkv_free_phages/quality_summary.tsv" | cut -f8 | sort | uniq -c | \
      awk '{printf "    %-20s: %s\n", $2, $1}'

    echo -n "  Total assessed: "
    tail -n +2 "$dir/phage_analysis/checkv_free_phages/quality_summary.tsv" | wc -l

    echo -n "  Mean completeness: "
    tail -n +2 "$dir/phage_analysis/checkv_free_phages/quality_summary.tsv" | \
      awk -F'\t' '{if ($9 != "NA") sum+=$9; count++} END {if (count>0) printf "%.1f%%\n", sum/count; else print "N/A"}'

    echo -n "  Mean contamination: "
    tail -n +2 "$dir/phage_analysis/checkv_free_phages/quality_summary.tsv" | \
      awk -F'\t' '{if ($10 != "NA") sum+=$10; count++} END {if (count>0) printf "%.1f%%\n", sum/count; else print "N/A"}'
  else
    echo "FREE PHAGES: CheckV results not found"
  fi
  echo ""
done

# Overall summary across all samples
echo "==============================================="
echo "OVERALL SUMMARY - ALL PROPHAGES (All Samples)"
echo "==============================================="

echo "Quality distribution:"
for dir in NovaSeq_*/; do
  if [ -f "$dir/phage_analysis/checkv_all_prophages/quality_summary.tsv" ]; then
    tail -n +2 "$dir/phage_analysis/checkv_all_prophages/quality_summary.tsv" | cut -f8
  fi
done | sort | uniq -c | awk '{printf "  %-20s: %s\n", $2, $1}'

echo -n "Total prophages assessed: "
for dir in NovaSeq_*/; do
  if [ -f "$dir/phage_analysis/checkv_all_prophages/quality_summary.tsv" ]; then
    tail -n +2 "$dir/phage_analysis/checkv_all_prophages/quality_summary.tsv"
  fi
done | wc -l

echo -n "Mean completeness: "
for dir in NovaSeq_*/; do
  if [ -f "$dir/phage_analysis/checkv_all_prophages/quality_summary.tsv" ]; then
    tail -n +2 "$dir/phage_analysis/checkv_all_prophages/quality_summary.tsv" | awk -F'\t' '{print $9}'
  fi
done | awk '{if ($1 != "NA") {sum+=$1; count++}} END {if (count>0) printf "%.1f%%\n", sum/count; else print "N/A"}'

echo -n "Mean contamination: "
for dir in NovaSeq_*/; do
  if [ -f "$dir/phage_analysis/checkv_all_prophages/quality_summary.tsv" ]; then
    tail -n +2 "$dir/phage_analysis/checkv_all_prophages/quality_summary.tsv" | awk -F'\t' '{print $10}'
  fi
done | awk '{if ($1 != "NA") {sum+=$1; count++}} END {if (count>0) printf "%.1f%%\n", sum/count; else print "N/A"}'

echo ""
echo "==============================================="
echo "OVERALL SUMMARY - FREE PHAGES (All Samples)"
echo "==============================================="

echo "Quality distribution:"
for dir in NovaSeq_*/; do
  if [ -f "$dir/phage_analysis/checkv_free_phages/quality_summary.tsv" ]; then
    tail -n +2 "$dir/phage_analysis/checkv_free_phages/quality_summary.tsv" | cut -f8
  fi
done | sort | uniq -c | awk '{printf "  %-20s: %s\n", $2, $1}'

echo -n "Total free phages assessed: "
for dir in NovaSeq_*/; do
  if [ -f "$dir/phage_analysis/checkv_free_phages/quality_summary.tsv" ]; then
    tail -n +2 "$dir/phage_analysis/checkv_free_phages/quality_summary.tsv"
  fi
done | wc -l

echo -n "Mean completeness: "
for dir in NovaSeq_*/; do
  if [ -f "$dir/phage_analysis/checkv_free_phages/quality_summary.tsv" ]; then
    tail -n +2 "$dir/phage_analysis/checkv_free_phages/quality_summary.tsv" | awk -F'\t' '{print $9}'
  fi
done | awk '{if ($1 != "NA") {sum+=$1; count++}} END {if (count>0) printf "%.1f%%\n", sum/count; else print "N/A"}'

echo -n "Mean contamination: "
for dir in NovaSeq_*/; do
  if [ -f "$dir/phage_analysis/checkv_free_phages/quality_summary.tsv" ]; then
    tail -n +2 "$dir/phage_analysis/checkv_free_phages/quality_summary.tsv" | awk -F'\t' '{print $10}'
  fi
done | awk '{if ($1 != "NA") {sum+=$1; count++}} END {if (count>0) printf "%.1f%%\n", sum/count; else print "N/A"}'
