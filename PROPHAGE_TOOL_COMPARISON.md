# Prophage Tool Comparison Pipeline

This branch (`prophage-tool-comparison`) extends the main pipeline to include PIDE as a third prophage detection tool and provides comprehensive comparison analysis between geNomad, PhiSpy, and PIDE.

## New Features

### 1. PIDE Integration
- **PIDE**: State-of-the-art prophage detection using ESM-2 protein language models
- Uses deep learning with 650 million parameters for superior boundary detection
- Excels at detecting longer prophages (>3kb) with precise coordinates

### 2. Comprehensive Tool Comparison
- Runs all three tools in parallel on the same input contigs
- Generates detailed comparison analysis without merging predictions
- Provides evidence-based data for integration decisions

## Usage

### Prerequisites
Before running the comparison pipeline, you'll need to:

1. **Download PIDE model**:
```bash
wget https://zenodo.org/records/12759619/files/PIDE.model.tar.gz
tar xzvf PIDE.model.tar.gz
```

2. **Clone PIDE repository**:
```bash
git clone https://github.com/chyghy/PIDE.git
```

3. **Update config.yaml paths**:
```yaml
pide_model: "/path/to/your/PIDE.model"     # Update with actual path
pide_script: "/path/to/your/PIDE"          # Update with actual path to PIDE repo
```

### Running the Comparison

#### Option 1: Run comparison analysis only
```bash
snakemake --profile ../profile/slurm/ --config reads=/path/to/reads outdir=/path/to/output run_prophage_comparison
```

#### Option 2: Run standard pipeline + comparison
```bash
# Run standard pipeline first
snakemake --profile ../profile/slurm/ --config reads=/path/to/reads outdir=/path/to/output run_everything

# Then run comparison
snakemake --profile ../profile/slurm/ --config reads=/path/to/reads outdir=/path/to/output run_prophage_comparison
```

## Output Files

The comparison analysis creates a new `comparison/` directory within each sample's `phage_analysis/` folder:

```
{sample}/phage_analysis/comparison/
├── raw_predictions.tsv       # All predictions from all tools with metadata
├── tool_statistics.tsv       # Summary stats per tool (counts, lengths, etc.)
├── overlap_analysis.tsv      # Detailed overlap analysis between tool pairs
├── agreement_summary.tsv     # Summary of tool agreement patterns
├── unique_predictions.tsv    # Predictions unique to each tool
├── unique_summary.tsv        # Summary of unique prediction counts
└── comparison_plots.pdf      # Visualization plots
```

## Key Analysis Outputs

### 1. Tool Statistics (`tool_statistics.tsv`)
- Total predictions per tool
- Length distributions (mean, median, min, max)
- Total base pairs predicted
- Number of unique contigs with predictions

### 2. Overlap Analysis (`overlap_analysis.tsv`)
- Pairwise overlaps between all tool predictions
- Overlap percentages for shared regions
- Boundary comparison data

### 3. Unique Predictions (`unique_predictions.tsv`)
- Prophages found by only one tool
- Tool-specific discovery patterns
- Potential complementary regions

### 4. Visualizations (`comparison_plots.pdf`)
- Prophage length distribution by tool
- Prediction count comparisons
- Overlap pattern visualizations

## Decision Framework

Use the comparison results to determine optimal integration strategy:

### High Complementarity Scenario
If tools find mostly unique prophages → Use union approach (combine all predictions)

### Tool Hierarchy Scenario  
If one tool consistently outperforms → Use hierarchical priority (e.g., geNomad > PIDE > PhiSpy)

### Length-Specific Scenario
If tools excel at different length ranges → Use hybrid strategy based on prophage length

## Example Analysis Questions

1. **Sensitivity**: Which tool finds the most prophages?
2. **Specificity**: Which predictions look most reliable?
3. **Complementarity**: Do tools find different prophage populations?
4. **Boundary accuracy**: Which tool provides better start/stop coordinates?
5. **Length bias**: Do tools have systematic length preferences?

## Technical Notes

### PIDE Requirements
- Python 3.8+
- PyTorch (CPU or GPU)
- fair-esm package
- 650M parameter ESM-2 model (~2.5GB download)

### Computational Considerations
- PIDE requires more memory than geNomad/PhiSpy
- ESM-2 model loading adds startup time
- Consider GPU acceleration for large datasets

### Known Limitations
- PIDE output format may need adjustment based on actual tool output
- Comparison assumes standardized contig naming (NODE_X format)
- Overlap tolerance currently set to exact coordinates (can be adjusted)

## Integration Recommendations

Based on comparison results, you can:
1. Update the main pipeline's merge strategy
2. Add PIDE as a third tool with appropriate priority
3. Implement length-based tool selection
4. Create consensus-based integration rules

## Next Steps

1. Run comparison on representative samples
2. Analyze results to determine optimal integration strategy
3. Update production pipeline based on findings
4. Consider merging improvements back to main branch