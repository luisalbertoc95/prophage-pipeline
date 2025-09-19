# Prophage Masking for MMseqs Taxonomy Implementation Plan

## Overview
Add optional prophage masking before MMseqs taxonomy to improve bacterial host classification by removing viral gene interference.

## Scientific Rationale
- Prophage regions contain viral genes that confuse MMseqs bacterial taxonomy
- Many contigs classified as "unknown/root" or "Caudoviricetes sp." could resolve to proper bacterial taxa
- Particularly valuable for unbinned contigs where GTDB-Tk isn't available

## Implementation Strategy

### 1. Configuration Parameter
Add to `config/config.yaml`:
```yaml
# Prophage masking for MMseqs taxonomy (optional)
mask_prophages_for_taxonomy: false  # Set to true to enable prophage masking before MMseqs
```

### 2. Command Line Usage
```bash
# Original approach (default)
snakemake --config mask_prophages_for_taxonomy=false

# With prophage masking
snakemake --config mask_prophages_for_taxonomy=true
```

### 3. Pipeline Modifications

#### A. Add Optional Masking Rule
```python
rule mask_prophages_for_mmseqs:
    input:
        contigs = os.path.join(config["outdir"], "{sample}", "binning", "final_filtered_contigs.fasta"),
        prophage_bed = os.path.join(config["outdir"], "{sample}", "phage_analysis", "merged_prophages.bed")
    output:
        masked_contigs = os.path.join(config["outdir"], "{sample}", "taxonomy", "contigs_masked_for_mmseqs.fasta")
    conda: config["conda_envs"]["bedtools"]  # Ensure bedtools is available
    shell:
        """
        bedtools subtract -a {input.contigs} -b {input.prophage_bed} > {output.masked_contigs}
        """
```

#### B. Conditional Input Function
```python
def get_mmseqs_input_contigs(wildcards):
    if config.get("mask_prophages_for_taxonomy", False):
        return os.path.join(config["outdir"], wildcards.sample, "taxonomy", "contigs_masked_for_mmseqs.fasta")
    else:
        return os.path.join(config["outdir"], wildcards.sample, "binning", "final_filtered_contigs.fasta")
```

#### C. Modify MMseqs Rule
Update the MMseqs taxonomy rule to use conditional input:
```python
rule mmseqs_taxonomy:
    input:
        contigs = get_mmseqs_input_contigs,
        db = config["mmseqs_database"]
    # ... rest of rule unchanged
```

#### D. Add Dependency Logic
Ensure masking rule runs when needed:
```python
def get_mmseqs_dependencies(wildcards):
    deps = []
    if config.get("mask_prophages_for_taxonomy", False):
        deps.append(os.path.join(config["outdir"], wildcards.sample, "taxonomy", "contigs_masked_for_mmseqs.fasta"))
    return deps
```

### 4. Technical Considerations

#### Fragment Handling
- `bedtools subtract` creates separate FASTA entries for bacterial regions flanking prophages
- Need to handle fragment-to-original contig mapping in downstream processing
- Multiple bacterial fragments per original contig may get different taxonomies

#### R Script Modifications
Current `add_taxonomy_to_prophage_table.R` expects one taxonomy per contig. May need to:
- Handle multiple MMseqs hits per original contig
- Take highest confidence assignment
- Or implement consensus logic

#### Minimum Length Filtering
- Some bacterial fragments may become too short for reliable taxonomy
- Consider adding minimum length filter after masking

### 5. Testing Strategy

#### Comparison Approach
1. **Baseline run**: `mask_prophages_for_taxonomy=false`
2. **Masked run**: `mask_prophages_for_taxonomy=true`
3. **Compare outputs**: Focus on reduction in "unknown/root" classifications

#### Key Metrics to Compare
- Number of "unknown/root" classifications
- Number of viral classifications (should decrease)
- Confidence scores for bacterial assignments
- Consistency with GTDB-Tk results for binned contigs

#### Test Samples
- Use samples with high prophage content (like `I14085_FGT_Metagenomic_FRESH_41680`)
- Focus on samples with many "unknown" classifications

### 6. Expected Outcomes
- Reduction in "unknown/root" classifications for unbinned contigs
- Fewer spurious viral classifications
- Higher confidence bacterial taxonomy assignments
- Better resolution of BV-associated bacteria (Prevotella, Peptoniphilus, etc.)

### 7. Implementation Priority
- **Phase 1**: Ensure current pipeline is robust and well-tested
- **Phase 2**: Implement optional prophage masking
- **Phase 3**: Comparative testing and validation

## Notes
- This is an optional enhancement - default behavior remains unchanged
- Adds minimal complexity while potentially providing significant scientific value
- Can be tested without disrupting current workflow