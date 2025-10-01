# Testing Checklist for Refactored Prophage Pipeline

This checklist verifies that the split MAG/unbinned workflow is functioning correctly.

## 1. Checkpoint & Per-MAG Processing

Check that bins were discovered and processed individually:

```bash
# Check bin discovery
cat out/{sample}/phage_analysis/mags/bin_list.txt

# Check per-MAG outputs exist (should have one directory per bin)
ls out/{sample}/phage_analysis/mags/bakta/
ls out/{sample}/phage_analysis/mags/phispy/
ls out/{sample}/phage_analysis/mags/genomad/
```

**What to verify:** Number of directories matches number of bins from DAS Tool

---

## 2. MAG Prophage Merging (bedtools logic)

```bash
# Check intermediate BED files
wc -l out/{sample}/phage_analysis/mags/genomad_all.bed
wc -l out/{sample}/phage_analysis/mags/phispy_all.bed
wc -l out/{sample}/phage_analysis/mags/phispy_unique.bed

# Check merged output
head out/{sample}/phage_analysis/mags/merged_prophages.bed
```

**What to verify:**
- `phispy_unique.bed` should have fewer lines than `phispy_all.bed` (overlaps removed)
- `merged_prophages.bed` = all genomad + unique phispy lines

---

## 3. Prophage Masking (unbinned)

```bash
# Check masking happened
cat out/{sample}/phage_analysis/unbinned/masking_stats.txt

# Look at masked contigs (should see N's where prophages were)
grep -A 5 "NODE_" out/{sample}/phage_analysis/unbinned/masked_contigs.fasta | head -20
```

**What to verify:**
- Masking stats shows number of regions masked
- Masked contigs contain stretches of N's

---

## 4. Dependency Order

Check log timestamps to verify execution order:

```bash
# MAG prophages should complete BEFORE gtdb-tk starts
ls -lht out/{sample}/phage_analysis/mags/prophage_table.tsv
ls -lht out/{sample}/taxonomy/gtdbtk/

# Masking should complete BEFORE mmseqs starts
ls -lht out/{sample}/phage_analysis/unbinned/masked_contigs.fasta
ls -lht out/{sample}/taxonomy/mmseqs/
```

**What to verify:** Prophage outputs have earlier timestamps than taxonomy

---

## 5. Final Output Structure

```bash
# Check all three FASTA outputs exist
ls -lh out/{sample}/phage_analysis/prophages_from_mags.fasta
ls -lh out/{sample}/phage_analysis/prophages_from_unbinned.fasta
ls -lh out/{sample}/phage_analysis/free_phages.fasta
ls -lh out/{sample}/phage_analysis/all_prophages_combined.fasta

# Check combined table has source column
head -1 out/{sample}/phage_analysis/final_prophage_table.tsv
# Should show: contig  start  end  tool  bin  source
```

**What to verify:**
- All three FASTAs exist and have content
- Table has 6 columns including "source"

---

## 6. Source Tracking in Tables

```bash
# Check MAG prophages are labeled correctly
head out/{sample}/phage_analysis/mags/prophage_table.tsv
# Source column should say "MAG", bin column should have bin.X

# Check unbinned prophages
head out/{sample}/phage_analysis/unbinned/prophage_table.tsv
# Source column should say "unbinned", bin column should have "none"

# Check merged table has both
awk '{print $6}' out/{sample}/phage_analysis/final_prophage_table.tsv | sort | uniq -c
# Should show counts for "source", "MAG", "unbinned"
```

---

## 7. Taxonomy Integration

```bash
# Check final table with taxonomy
head out/{sample}/phage_analysis/final_prophage_table_with_host_taxonomy.tsv

# Verify MAG prophages have GTDB-Tk taxonomy
grep "MAG" out/{sample}/phage_analysis/final_prophage_table_with_host_taxonomy.tsv | head -5

# Verify unbinned prophages have MMseqs taxonomy
grep "unbinned" out/{sample}/phage_analysis/final_prophage_table_with_host_taxonomy.tsv | head -5
```

**What to verify:** Both MAG and unbinned prophages have appropriate taxonomy annotations

---

## 8. CheckV Results

```bash
# Should include both prophages AND free phages
wc -l out/{sample}/phage_analysis/unbinned/checkv/quality_summary.tsv
```

---

## Quick Summary Commands

One-liner to check key outputs:

```bash
sample="your_sample_name"

echo "=== Bins processed ===" && cat out/${sample}/phage_analysis/mags/bin_list.txt | wc -l

echo "=== MAG prophages ===" && grep -c "^>" out/${sample}/phage_analysis/prophages_from_mags.fasta

echo "=== Unbinned prophages ===" && grep -c "^>" out/${sample}/phage_analysis/prophages_from_unbinned.fasta

echo "=== Free phages ===" && grep -c "^>" out/${sample}/phage_analysis/free_phages.fasta

echo "=== Source distribution ===" && tail -n +2 out/${sample}/phage_analysis/final_prophage_table.tsv | awk '{print $6}' | sort | uniq -c
```

---

## Key Success Criteria

✅ **Per-MAG processing worked**: Each bin has its own bakta/phispy/genomad output directory

✅ **Bedtools merging worked**: phispy_unique.bed has fewer entries than phispy_all.bed

✅ **Masking worked**: masked_contigs.fasta contains N's in prophage regions

✅ **Dependencies correct**: Prophage prediction completed before taxonomy assignment

✅ **Three output types**: Separate FASTAs for MAG prophages, unbinned prophages, and free phages

✅ **Source tracking**: Final table correctly labels MAG vs unbinned prophages

✅ **Taxonomy**: MAGs have GTDB-Tk taxonomy, unbinned contigs have MMseqs taxonomy
