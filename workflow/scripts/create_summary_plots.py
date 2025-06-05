#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
import argparse
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

def collect_prophage_data(outdir):
    """Collect prophage data from all samples"""
    print("Collecting prophage data...")
    
    # Find all prophage table files
    prophage_files = glob.glob(
        os.path.join(outdir, "phage_analysis", "*", "final_prophage_table_with_host_taxonomy.tsv")
    )
    
    if not prophage_files:
        raise ValueError("No prophage data found. Please check that the pipeline has completed successfully.")
    
    all_data = []
    for file_path in prophage_files:
        # Extract sample name from path
        sample = os.path.basename(os.path.dirname(file_path))
        
        try:
            df = pd.read_csv(file_path, sep='\t')
            df['sample'] = sample
            all_data.append(df)
        except Exception as e:
            print(f"Warning: Could not read {file_path}: {e}")
    
    if not all_data:
        raise ValueError("No valid prophage data could be read.")
    
    return pd.concat(all_data, ignore_index=True)

def collect_checkv_data(outdir):
    """Collect CheckV quality data"""
    print("Collecting CheckV quality data...")
    
    checkv_files = glob.glob(
        os.path.join(outdir, "phage_analysis", "*", "checkv", "quality_summary.tsv")
    )
    
    all_data = []
    for file_path in checkv_files:
        sample = os.path.basename(os.path.dirname(os.path.dirname(file_path)))
        
        try:
            df = pd.read_csv(file_path, sep='\t')
            df['sample'] = sample
            all_data.append(df)
        except Exception as e:
            print(f"Warning: Could not read {file_path}: {e}")
    
    return pd.concat(all_data, ignore_index=True) if all_data else pd.DataFrame()

def create_summary_plots(prophage_data, checkv_data, output_dir):
    """Create summary visualizations"""
    print("Creating visualizations...")
    
    # Set up plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create output directory
    plot_dir = os.path.join(output_dir, "summary_plots")
    os.makedirs(plot_dir, exist_ok=True)
    
    # 1. Prophages per sample
    plt.figure(figsize=(12, 6))
    sample_counts = prophage_data['sample'].value_counts()
    sample_counts.plot(kind='bar')
    plt.title('Total Prophages Detected per Sample')
    plt.xlabel('Sample')
    plt.ylabel('Number of Prophages')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'prophages_per_sample.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Tool comparison
    if 'tool' in prophage_data.columns:
        plt.figure(figsize=(8, 6))
        tool_counts = prophage_data['tool'].value_counts()
        colors = ['steelblue', 'orange']
        tool_counts.plot(kind='bar', color=colors[:len(tool_counts)])
        plt.title('Prophage Detection by Tool')
        plt.xlabel('Detection Tool')
        plt.ylabel('Number of Prophages')
        plt.xticks(rotation=45, ha='right')
        
        # Add percentage labels
        total = tool_counts.sum()
        for i, v in enumerate(tool_counts.values):
            plt.text(i, v + total*0.01, f'{v}\n({v/total*100:.1f}%)', 
                    ha='center', va='bottom')
        
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'tool_comparison.png'), dpi=300, bbox_inches='tight')
        plt.close()
    
    # 3. Tool comparison by sample (stacked bar)
    if 'tool' in prophage_data.columns:
        plt.figure(figsize=(12, 6))
        tool_by_sample = prophage_data.groupby(['sample', 'tool']).size().unstack(fill_value=0)
        tool_by_sample.plot(kind='bar', stacked=True, figsize=(12, 6))
        plt.title('Prophage Detection Tools by Sample')
        plt.xlabel('Sample')
        plt.ylabel('Number of Prophages')
        plt.xticks(rotation=45, ha='right')
        plt.legend(title='Tool')
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'tool_by_sample.png'), dpi=300, bbox_inches='tight')
        plt.close()
    
    # 4. Host taxonomy distribution (top 10)
    if 'phylum' in prophage_data.columns:
        plt.figure(figsize=(10, 6))
        phylum_counts = prophage_data['phylum'].value_counts().head(10)
        phylum_counts.plot(kind='barh', color='forestgreen')
        plt.title('Top 10 Host Phyla for Prophages')
        plt.xlabel('Number of Prophages')
        plt.ylabel('Phylum')
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'taxonomy_distribution.png'), dpi=300, bbox_inches='tight')
        plt.close()
    
    # 5. CheckV quality distribution
    if not checkv_data.empty and 'checkv_quality' in checkv_data.columns:
        plt.figure(figsize=(10, 6))
        quality_counts = checkv_data['checkv_quality'].value_counts()
        quality_counts.plot(kind='bar', color='purple', alpha=0.7)
        plt.title('Prophage Quality Distribution (CheckV)')
        plt.xlabel('Quality Category')
        plt.ylabel('Number of Prophages')
        plt.xticks(rotation=45, ha='right')
        
        # Add percentage labels
        total = quality_counts.sum()
        for i, v in enumerate(quality_counts.values):
            plt.text(i, v + total*0.01, f'{v}\n({v/total*100:.1f}%)', 
                    ha='center', va='bottom')
        
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'quality_distribution.png'), dpi=300, bbox_inches='tight')
        plt.close()
    
    # 6. Prophage length distribution
    if not checkv_data.empty and 'provirus_length' in checkv_data.columns:
        plt.figure(figsize=(10, 6))
        plt.hist(checkv_data['provirus_length'], bins=30, color='orange', alpha=0.7, edgecolor='black')
        plt.title('Prophage Length Distribution')
        plt.xlabel('Prophage Length (bp)')
        plt.ylabel('Count')
        plt.xscale('log')
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'length_distribution.png'), dpi=300, bbox_inches='tight')
        plt.close()
    
    print(f"Plots saved to: {plot_dir}")
    return plot_dir

def create_summary_tables(prophage_data, checkv_data, output_dir):
    """Create summary tables"""
    print("Creating summary tables...")
    
    # Sample summary
    sample_summary = prophage_data.groupby('sample').agg({
        'contig': 'count',  # total prophages
        'tool': lambda x: (x == 'genomad').sum() if 'tool' in prophage_data.columns else 0,
        'contig': ['count', 'nunique']  # total and unique contigs
    }).round(2)
    
    # Flatten column names
    sample_summary.columns = ['total_prophages', 'genomad_prophages', 'total_prophages_2', 'unique_contigs']
    sample_summary = sample_summary[['total_prophages', 'genomad_prophages', 'unique_contigs']]
    sample_summary['phispy_prophages'] = sample_summary['total_prophages'] - sample_summary['genomad_prophages']
    
    # Tool summary
    if 'tool' in prophage_data.columns:
        tool_summary = prophage_data['tool'].value_counts().reset_index()
        tool_summary.columns = ['tool', 'prophages_detected']
        tool_summary['percentage'] = (tool_summary['prophages_detected'] / tool_summary['prophages_detected'].sum() * 100).round(1)
    else:
        tool_summary = pd.DataFrame(columns=['tool', 'prophages_detected', 'percentage'])
    
    # Save tables
    sample_summary.to_csv(os.path.join(output_dir, 'prophage_summary_by_sample.tsv'), sep='\t')
    tool_summary.to_csv(os.path.join(output_dir, 'tool_detection_summary.tsv'), sep='\t', index=False)
    
    print(f"Summary tables saved to: {output_dir}")
    return sample_summary, tool_summary

def main(output_directory):
    """Main analysis function"""
    print("Starting prophage pipeline summary analysis...")
    print(f"Output directory: {output_directory}")
    
    # Collect data
    prophage_data = collect_prophage_data(output_directory)
    checkv_data = collect_checkv_data(output_directory)
    
    # Create summaries
    sample_summary, tool_summary = create_summary_tables(prophage_data, checkv_data, output_directory)
    plot_dir = create_summary_plots(prophage_data, checkv_data, output_directory)
    
    # Print basic statistics
    print("\n=== SUMMARY STATISTICS ===")
    print(f"Total samples: {prophage_data['sample'].nunique()}")
    print(f"Total prophages detected: {len(prophage_data)}")
    print(f"Average prophages per sample: {len(prophage_data) / prophage_data['sample'].nunique():.1f}")
    
    if 'tool' in prophage_data.columns:
        print("\nDetection tool breakdown:")
        for _, row in tool_summary.iterrows():
            print(f"  {row['tool']}: {row['prophages_detected']} ({row['percentage']}%)")
    
    print(f"\nResults saved to: {output_directory}")
    print(f"Plots directory: {plot_dir}")
    
    return {
        'prophage_data': prophage_data,
        'checkv_data': checkv_data,
        'sample_summary': sample_summary,
        'tool_summary': tool_summary,
        'plot_dir': plot_dir
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create summary plots and statistics for prophage pipeline results')
    parser.add_argument('output_dir', help='Output directory containing pipeline results')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.output_dir):
        raise ValueError(f"Output directory does not exist: {args.output_dir}")
    
    main(args.output_dir)