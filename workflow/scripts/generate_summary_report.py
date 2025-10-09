#!/usr/bin/env python3
"""
Generate comprehensive HTML summary report for prophage pipeline results.
Aggregates data from all samples with statistical analysis and interactive visualizations.
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from Bio import SeqIO
from pathlib import Path
import sys
from collections import defaultdict

# ============================================================================
# STATISTICAL FUNCTIONS
# ============================================================================

def calculate_stats(values):
    """Calculate comprehensive statistics for a list of values"""
    if not values or len(values) == 0:
        return {
            'total': 0, 'mean': 0, 'median': 0, 'std': 0,
            'min': 0, 'max': 0, 'range': 0, 'count': 0
        }

    return {
        'total': sum(values),
        'mean': np.mean(values),
        'median': np.median(values),
        'std': np.std(values),
        'min': np.min(values),
        'max': np.max(values),
        'range': np.max(values) - np.min(values),
        'count': len(values)
    }

def format_stat(value, decimals=1):
    """Format a statistic for display"""
    if isinstance(value, (int, np.integer)):
        return f"{value:,}"
    else:
        return f"{value:,.{decimals}f}"

# ============================================================================
# DATA PARSING FUNCTIONS
# ============================================================================

def parse_prophage_tables(samples, outdir):
    """Parse all prophage tables and aggregate statistics"""
    data = []

    for sample in samples:
        sample_data = {'sample': sample}

        # Basic prophage table (has all prophages with source info)
        basic_table_path = Path(outdir) / sample / "phage_analysis" / "final_prophage_table.tsv"
        if basic_table_path.exists():
            df = pd.read_csv(basic_table_path, sep='\t')
            sample_data['total_prophages'] = len(df)
            sample_data['mag_prophages'] = len(df[df['source'] == 'MAG'])
            sample_data['unbinned_prophages'] = len(df[df['source'] == 'unbinned'])
            sample_data['genomad_prophages'] = len(df[df['tool'] == 'genomad'])
            sample_data['phispy_prophages'] = len(df[df['tool'] == 'phispy'])
        else:
            sample_data.update({
                'total_prophages': 0, 'mag_prophages': 0, 'unbinned_prophages': 0,
                'genomad_prophages': 0, 'phispy_prophages': 0
            })

        # MAG table for bin info
        mag_table_path = Path(outdir) / sample / "phage_analysis" / "mags" / "prophage_table.tsv"
        if mag_table_path.exists():
            mag_df = pd.read_csv(mag_table_path, sep='\t')
            sample_data['mags_with_prophages'] = mag_df['bin'].nunique() if len(mag_df) > 0 else 0
        else:
            sample_data['mags_with_prophages'] = 0

        # Taxonomy table (when available)
        tax_table_path = Path(outdir) / sample / "phage_analysis" / "final_prophage_table_with_host_taxonomy.tsv"
        if tax_table_path.exists():
            tax_df = pd.read_csv(tax_table_path, sep='\t')
            sample_data['prophages_with_taxonomy'] = len(tax_df)
        else:
            sample_data['prophages_with_taxonomy'] = 0

        data.append(sample_data)

    return pd.DataFrame(data)

def parse_free_phages(samples, outdir):
    """Parse free phage detection results"""
    data = []

    for sample in samples:
        free_phage_table = Path(outdir) / sample / "phage_analysis" / "unbinned" / "free_phage_table.tsv"
        if free_phage_table.exists():
            df = pd.read_csv(free_phage_table, sep='\t')
            free_phage_count = len(df)
        else:
            free_phage_count = 0

        data.append({'sample': sample, 'free_phages': free_phage_count})

    return pd.DataFrame(data)

def parse_checkv_quality(samples, outdir, checkv_type='all_prophages'):
    """Parse CheckV quality assessment results"""
    data = []

    for sample in samples:
        checkv_path = Path(outdir) / sample / "phage_analysis" / f"checkv_{checkv_type}" / "quality_summary.tsv"
        if checkv_path.exists() and checkv_path.stat().st_size > 0:
            try:
                df = pd.read_csv(checkv_path, sep='\t')

                for _, row in df.iterrows():
                    data.append({
                        'sample': sample,
                        'contig': row.get('contig_id', 'unknown'),
                        'checkv_quality': row.get('checkv_quality', 'Not-determined'),
                        'completeness': row.get('completeness', 0),
                        'contamination': row.get('contamination', 0)
                    })
            except Exception as e:
                print(f"Warning: Could not parse CheckV for {sample}: {e}", file=sys.stderr)

    return pd.DataFrame(data) if data else pd.DataFrame(columns=['sample', 'contig', 'checkv_quality', 'completeness', 'contamination'])

def parse_checkm_quality(samples, outdir):
    """Parse CheckM bin quality results (excludes unbinned contigs)"""
    data = []

    for sample in samples:
        checkm_path = Path(outdir) / sample / "binning" / "checkm" / "checkm_out.tsv"
        if checkm_path.exists() and checkm_path.stat().st_size > 0:
            try:
                df = pd.read_csv(checkm_path, sep='\t')

                for _, row in df.iterrows():
                    # CheckM output columns may vary, adapt as needed
                    bin_id = row.get('Bin Id', row.iloc[0] if len(row) > 0 else 'unknown')

                    # Skip unbinned contigs (identified by NODE_ in name)
                    if 'NODE_' in str(bin_id):
                        continue

                    completeness = row.get('Completeness', row.iloc[13] if len(row) > 13 else 0)
                    contamination = row.get('Contamination', row.iloc[14] if len(row) > 14 else 0)

                    data.append({
                        'sample': sample,
                        'bin': bin_id,
                        'completeness': float(completeness),
                        'contamination': float(contamination)
                    })
            except Exception as e:
                print(f"Warning: Could not parse CheckM for {sample}: {e}", file=sys.stderr)

    return pd.DataFrame(data) if data else pd.DataFrame(columns=['sample', 'bin', 'completeness', 'contamination'])

def calculate_assembly_stats(fasta_path):
    """Calculate assembly statistics from FASTA file"""
    try:
        sequences = list(SeqIO.parse(fasta_path, "fasta"))
        if not sequences:
            return {'n_contigs': 0, 'n50': 0, 'longest_contig': 0, 'total_length': 0}

        lengths = sorted([len(seq) for seq in sequences], reverse=True)
        total_length = sum(lengths)

        # Calculate N50
        cumsum = 0
        n50 = 0
        for length in lengths:
            cumsum += length
            if cumsum >= total_length / 2:
                n50 = length
                break

        return {
            'n_contigs': len(sequences),
            'n50': n50,
            'longest_contig': lengths[0],
            'total_length': total_length
        }
    except Exception as e:
        print(f"Warning: Could not parse assembly {fasta_path}: {e}", file=sys.stderr)
        return {'n_contigs': 0, 'n50': 0, 'longest_contig': 0, 'total_length': 0}

def parse_assembly_stats(samples, outdir):
    """Parse assembly statistics for all samples"""
    data = []

    for sample in samples:
        assembly_path = Path(outdir) / sample / "assembly" / "final.contigs.fa"
        stats = calculate_assembly_stats(assembly_path)
        stats['sample'] = sample
        data.append(stats)

    return pd.DataFrame(data)

def parse_taxonomy(samples, outdir):
    """Parse taxonomy assignments from prophage tables"""
    all_taxonomy = []

    for sample in samples:
        tax_path = Path(outdir) / sample / "phage_analysis" / "final_prophage_table_with_host_taxonomy.tsv"
        if tax_path.exists():
            try:
                df = pd.read_csv(tax_path, sep='\t')
                if 'phylum' in df.columns and 'class' in df.columns:
                    for _, row in df.iterrows():
                        all_taxonomy.append({
                            'sample': sample,
                            'phylum': row.get('phylum', 'Unknown'),
                            'class': row.get('class', 'Unknown'),
                            'taxonomy_source': row.get('taxonomy_source', 'Unknown')
                        })
            except Exception as e:
                print(f"Warning: Could not parse taxonomy for {sample}: {e}", file=sys.stderr)

    return pd.DataFrame(all_taxonomy) if all_taxonomy else pd.DataFrame(columns=['sample', 'phylum', 'class', 'taxonomy_source'])

# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

def create_summary_cards(prophage_df, free_phage_df, checkm_df, assembly_df):
    """Create HTML summary statistics cards"""

    # Calculate statistics
    prophage_stats = calculate_stats(prophage_df['total_prophages'].tolist())
    mag_prophage_stats = calculate_stats(prophage_df['mag_prophages'].tolist())
    unbinned_prophage_stats = calculate_stats(prophage_df['unbinned_prophages'].tolist())
    free_phage_stats = calculate_stats(free_phage_df['free_phages'].tolist())
    mags_stats = calculate_stats(checkm_df.groupby('sample').size().tolist()) if len(checkm_df) > 0 else calculate_stats([])
    mags_with_prophages_stats = calculate_stats(prophage_df['mags_with_prophages'].tolist())
    n50_stats = calculate_stats(assembly_df['n50'].tolist())

    # High-quality MAGs (>90% complete, <5% contamination)
    if len(checkm_df) > 0:
        hq_mags = len(checkm_df[(checkm_df['completeness'] >= 90) & (checkm_df['contamination'] < 5)])
    else:
        hq_mags = 0

    cards_html = f"""
    <div class="summary-cards">
        <div class="card">
            <h3>Total Samples</h3>
            <div class="stat-value">{len(prophage_df)}</div>
        </div>
        <div class="card">
            <h3>Total Prophages</h3>
            <div class="stat-value">{format_stat(prophage_stats['total'])}</div>
            <div class="stat-detail">Across all samples</div>
        </div>
        <div class="card">
            <h3>MAG Prophages</h3>
            <div class="stat-value">{format_stat(mag_prophage_stats['total'])}</div>
            <div class="stat-detail">Across all samples</div>
        </div>
        <div class="card">
            <h3>Unbinned Prophages</h3>
            <div class="stat-value">{format_stat(unbinned_prophage_stats['total'])}</div>
            <div class="stat-detail">Across all samples</div>
        </div>
        <div class="card">
            <h3>Free Phages</h3>
            <div class="stat-value">{format_stat(free_phage_stats['total'])}</div>
            <div class="stat-detail">Across all samples</div>
        </div>
        <div class="card">
            <h3>MAGs with Prophages</h3>
            <div class="stat-value">{format_stat(mags_with_prophages_stats['total'])}</div>
            <div class="stat-detail">Total MAGs: {format_stat(len(checkm_df))}</div>
        </div>
        <div class="card">
            <h3>High-Quality MAGs</h3>
            <div class="stat-value">{hq_mags}</div>
            <div class="stat-detail">>90% complete, <5% contamination</div>
        </div>
    </div>
    """

    return cards_html

def create_stats_table(prophage_df, free_phage_df, checkv_all_df, checkv_free_df, checkm_df, assembly_df):
    """Create comprehensive statistical summary table"""

    rows = []

    # Prophage metrics
    rows.append(create_stats_row("Prophages", prophage_df['total_prophages'].tolist()))
    rows.append(create_stats_row("MAG prophages", prophage_df['mag_prophages'].tolist()))
    rows.append(create_stats_row("Unbinned prophages", prophage_df['unbinned_prophages'].tolist()))
    rows.append(create_stats_row("Free phages", free_phage_df['free_phages'].tolist()))

    # Assembly metrics
    rows.append(create_stats_row("Total contigs", assembly_df['n_contigs'].tolist()))
    rows.append(create_stats_row("Assembly N50 (bp)", assembly_df['n50'].tolist()))
    rows.append(create_stats_row("Total assembly length (bp)", assembly_df['total_length'].tolist()))

    # MAG metrics
    if len(checkm_df) > 0:
        mags_per_sample = checkm_df.groupby('sample').size().tolist()
        rows.append(create_stats_row("MAGs", mags_per_sample))
        rows.append(create_stats_row("MAG completeness (%) [CheckM]", checkm_df['completeness'].tolist(), decimals=1))
        rows.append(create_stats_row("MAG contamination (%) [CheckM]", checkm_df['contamination'].tolist(), decimals=1))

    # CheckV metrics for prophages
    if len(checkv_all_df) > 0:
        checkv_comp = checkv_all_df[checkv_all_df['completeness'] > 0]['completeness'].tolist()
        checkv_cont = checkv_all_df[checkv_all_df['contamination'] > 0]['contamination'].tolist()
        if checkv_comp:
            rows.append(create_stats_row("Prophage completeness (%) [CheckV]", checkv_comp, decimals=1))
        if checkv_cont:
            rows.append(create_stats_row("Prophage contamination (%) [CheckV]", checkv_cont, decimals=1))

    # CheckV metrics for free phages
    if len(checkv_free_df) > 0:
        free_checkv_comp = checkv_free_df[checkv_free_df['completeness'] > 0]['completeness'].tolist()
        free_checkv_cont = checkv_free_df[checkv_free_df['contamination'] > 0]['contamination'].tolist()
        if free_checkv_comp:
            rows.append(create_stats_row("Free phage completeness (%) [CheckV]", free_checkv_comp, decimals=1))
        if free_checkv_cont:
            rows.append(create_stats_row("Free phage contamination (%) [CheckV]", free_checkv_cont, decimals=1))

    table_html = """
    <table class="stats-table">
        <thead>
            <tr>
                <th>Metric (per sample)</th>
                <th>Mean</th>
                <th>Median</th>
                <th>StdDev</th>
                <th>Min</th>
                <th>Max</th>
                <th>Range</th>
            </tr>
        </thead>
        <tbody>
    """

    table_html += "\n".join(rows)
    table_html += """
        </tbody>
    </table>
    """

    return table_html

def create_stats_row(metric_name, values, decimals=0):
    """Create a single row for the statistics table

    Args:
        metric_name: Name of the metric
        values: List of values (per sample)
        decimals: Number of decimal places for formatting
    """
    stats = calculate_stats(values)

    return f"""
            <tr>
                <td>{metric_name}</td>
                <td>{format_stat(stats['mean'], decimals)}</td>
                <td>{format_stat(stats['median'], decimals)}</td>
                <td>{format_stat(stats['std'], decimals)}</td>
                <td>{format_stat(stats['min'], decimals)}</td>
                <td>{format_stat(stats['max'], decimals)}</td>
                <td>{format_stat(stats['range'], decimals)}</td>
            </tr>
    """

def plot_prophages_per_sample(prophage_df):
    """Create stacked bar chart of prophages per sample"""

    fig = go.Figure(data=[
        go.Bar(name='MAG', x=prophage_df['sample'], y=prophage_df['mag_prophages']),
        go.Bar(name='Unbinned', x=prophage_df['sample'], y=prophage_df['unbinned_prophages'])
    ])

    fig.update_layout(
        title='Prophages Per Sample (MAG vs Unbinned)',
        barmode='stack',
        xaxis_title='Sample',
        yaxis_title='Number of Prophages',
        height=500
    )

    return fig.to_html(full_html=False, include_plotlyjs='cdn')

def plot_prophage_box(prophage_df):
    """Create box plot of prophage distribution"""

    fig = go.Figure()

    fig.add_trace(go.Box(
        y=prophage_df['total_prophages'],
        x=['Total Prophages'] * len(prophage_df),
        name='Total'
    ))

    fig.add_trace(go.Box(
        y=prophage_df['mag_prophages'],
        x=['MAG Prophages'] * len(prophage_df),
        name='MAG'
    ))

    fig.add_trace(go.Box(
        y=prophage_df['unbinned_prophages'],
        x=['Unbinned Prophages'] * len(prophage_df),
        name='Unbinned'
    ))

    fig.update_layout(
        title='Prophage Distribution Across Samples',
        yaxis_title='Number of Prophages',
        height=500,
        showlegend=False
    )

    return fig.to_html(full_html=False, include_plotlyjs='cdn')

def plot_detection_tools(prophage_df):
    """Create stacked bar chart of detection tool breakdown"""

    fig = go.Figure(data=[
        go.Bar(name='GeNomad', x=prophage_df['sample'], y=prophage_df['genomad_prophages']),
        go.Bar(name='PhiSpy (unique)', x=prophage_df['sample'], y=prophage_df['phispy_prophages'])
    ])

    fig.update_layout(
        title='Prophage Detection Tool Breakdown',
        barmode='stack',
        xaxis_title='Sample',
        yaxis_title='Number of Prophages',
        height=500
    )

    return fig.to_html(full_html=False, include_plotlyjs='cdn')

def plot_checkv_quality(checkv_df):
    """Create stacked bar chart of CheckV quality tiers"""

    if len(checkv_df) == 0:
        return "<p>No CheckV data available</p>"

    # Count quality tiers per sample
    quality_counts = checkv_df.groupby(['sample', 'checkv_quality']).size().unstack(fill_value=0)

    # Create stacked bar chart
    fig = go.Figure()

    quality_order = ['Complete', 'High-quality', 'Medium-quality', 'Low-quality', 'Not-determined']
    colors = ['#2ecc71', '#3498db', '#f39c12', '#e74c3c', '#95a5a6']

    for quality, color in zip(quality_order, colors):
        if quality in quality_counts.columns:
            fig.add_trace(go.Bar(
                name=quality,
                x=quality_counts.index,
                y=quality_counts[quality],
                marker_color=color
            ))

    fig.update_layout(
        title='CheckV Quality Assessment (All Prophages)',
        barmode='stack',
        xaxis_title='Sample',
        yaxis_title='Number of Prophages',
        height=500
    )

    return fig.to_html(full_html=False, include_plotlyjs='cdn')

def plot_checkv_completeness_box(checkv_df):
    """Create box plot of CheckV completeness"""

    if len(checkv_df) == 0 or checkv_df['completeness'].sum() == 0:
        return "<p>No CheckV completeness data available</p>"

    # Filter out zero completeness
    data = checkv_df[checkv_df['completeness'] > 0]

    fig = go.Figure(data=[go.Box(
        y=data['completeness'],
        name='Completeness',
        marker_color='#3498db'
    )])

    fig.update_layout(
        title='Prophage Completeness Distribution (CheckV)',
        yaxis_title='Completeness (%)',
        height=400
    )

    return fig.to_html(full_html=False, include_plotlyjs='cdn')

def plot_checkm_scatter(checkm_df):
    """Create scatter plot of CheckM completeness vs contamination"""

    if len(checkm_df) == 0:
        return "<p>No CheckM data available</p>"

    # Determine quality tier
    def get_quality_tier(row):
        if row['completeness'] >= 90 and row['contamination'] < 5:
            return 'High-quality (>90% complete, <5% contam)'
        elif row['completeness'] >= 50 and row['contamination'] < 10:
            return 'Medium-quality (>50% complete, <10% contam)'
        else:
            return 'Low-quality'

    checkm_df['quality_tier'] = checkm_df.apply(get_quality_tier, axis=1)

    fig = px.scatter(
        checkm_df,
        x='completeness',
        y='contamination',
        color='quality_tier',
        hover_data=['sample', 'bin'],
        title='MAG Quality (CheckM)',
        labels={'completeness': 'Completeness (%)', 'contamination': 'Contamination (%)'},
        color_discrete_map={
            'High-quality (>90% complete, <5% contam)': '#2ecc71',
            'Medium-quality (>50% complete, <10% contam)': '#f39c12',
            'Low-quality': '#e74c3c'
        }
    )

    fig.update_layout(height=500)

    return fig.to_html(full_html=False, include_plotlyjs='cdn')

def plot_checkm_box(checkm_df):
    """Create box plots for CheckM completeness and contamination"""

    if len(checkm_df) == 0:
        return "<p>No CheckM data available</p>"

    fig = make_subplots(rows=1, cols=2, subplot_titles=('Completeness', 'Contamination'))

    fig.add_trace(
        go.Box(y=checkm_df['completeness'], name='Completeness', marker_color='#3498db'),
        row=1, col=1
    )

    fig.add_trace(
        go.Box(y=checkm_df['contamination'], name='Contamination', marker_color='#e74c3c'),
        row=1, col=2
    )

    fig.update_xaxes(title_text="", row=1, col=1)
    fig.update_xaxes(title_text="", row=1, col=2)
    fig.update_yaxes(title_text="Percentage (%)", row=1, col=1)

    fig.update_layout(
        title_text='MAG Quality Distribution (CheckM)',
        height=400,
        showlegend=False
    )

    return fig.to_html(full_html=False, include_plotlyjs='cdn')

def plot_taxonomy_phylum(taxonomy_df):
    """Create bar chart of host taxonomy at phylum level"""

    if len(taxonomy_df) == 0:
        return "<p>No taxonomy data available</p>"

    # Count phyla
    phylum_counts = taxonomy_df['phylum'].value_counts().head(10)

    fig = go.Figure(data=[
        go.Bar(x=phylum_counts.index, y=phylum_counts.values, marker_color='#9b59b6')
    ])

    fig.update_layout(
        title='Top 10 Prophage Host Taxa (Phylum Level)',
        xaxis_title='Phylum',
        yaxis_title='Number of Prophages',
        height=500
    )

    return fig.to_html(full_html=False, include_plotlyjs='cdn')

def plot_taxonomy_source(taxonomy_df):
    """Create pie chart of taxonomy source (GTDB-Tk vs MMseqs)"""

    if len(taxonomy_df) == 0:
        return "<p>No taxonomy data available</p>"

    source_counts = taxonomy_df['taxonomy_source'].value_counts()

    fig = go.Figure(data=[go.Pie(
        labels=source_counts.index,
        values=source_counts.values,
        hole=0.3
    )])

    fig.update_layout(
        title='Taxonomy Assignment Source',
        height=400
    )

    return fig.to_html(full_html=False, include_plotlyjs='cdn')

def plot_assembly_n50_box(assembly_df):
    """Create box plot of assembly N50"""

    fig = go.Figure(data=[go.Box(
        y=assembly_df['n50'] / 1000,  # Convert to kb
        name='N50',
        marker_color='#1abc9c'
    )])

    fig.update_layout(
        title='Assembly N50 Distribution',
        yaxis_title='N50 (kb)',
        height=400
    )

    return fig.to_html(full_html=False, include_plotlyjs='cdn')

# ============================================================================
# HTML GENERATION
# ============================================================================

def generate_html_report(prophage_df, free_phage_df, checkv_all_df, checkv_free_df, checkm_df, assembly_df, taxonomy_df, output_path):
    """Generate complete HTML report"""

    # Generate all plots
    summary_cards = create_summary_cards(prophage_df, free_phage_df, checkm_df, assembly_df)
    stats_table = create_stats_table(prophage_df, free_phage_df, checkv_all_df, checkv_free_df, checkm_df, assembly_df)

    prophage_bar = plot_prophages_per_sample(prophage_df)
    prophage_box = plot_prophage_box(prophage_df)
    detection_tools = plot_detection_tools(prophage_df)

    checkv_quality = plot_checkv_quality(checkv_all_df)
    checkv_comp_box = plot_checkv_completeness_box(checkv_all_df)

    checkm_scatter = plot_checkm_scatter(checkm_df)
    checkm_box = plot_checkm_box(checkm_df)

    taxonomy_phylum = plot_taxonomy_phylum(taxonomy_df)
    taxonomy_source = plot_taxonomy_source(taxonomy_df)

    n50_box = plot_assembly_n50_box(assembly_df)

    # CSS styling
    css = """
    <style>
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }
        .container {
            max-width: 1400px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        h1 {
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }
        h2 {
            color: #34495e;
            margin-top: 40px;
            border-left: 4px solid #3498db;
            padding-left: 15px;
        }
        .summary-cards {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 30px 0;
        }
        .card {
            background-color: #667eea;
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }
        .card h3 {
            margin: 0 0 10px 0;
            font-size: 14px;
            opacity: 0.9;
        }
        .stat-value {
            font-size: 32px;
            font-weight: bold;
            margin: 10px 0;
        }
        .stat-detail {
            font-size: 12px;
            opacity: 0.8;
            margin: 5px 0;
        }
        .stats-table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            font-size: 14px;
        }
        .stats-table th {
            background-color: #3498db;
            color: white;
            padding: 12px;
            text-align: left;
        }
        .stats-table td {
            padding: 10px;
            border-bottom: 1px solid #ddd;
        }
        .stats-table tr:hover {
            background-color: #f5f5f5;
        }
        .plot-container {
            margin: 30px 0;
            padding: 20px;
            background-color: #fafafa;
            border-radius: 8px;
        }
        .timestamp {
            color: #7f8c8d;
            font-size: 12px;
            margin-top: 40px;
            text-align: center;
        }
    </style>
    """

    # Build HTML
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="UTF-8">
        <title>Prophage Pipeline Summary Report</title>
        {css}
    </head>
    <body>
        <div class="container">
            <h1>Prophage Pipeline Summary Report</h1>

            <h2>Overall Totals</h2>
            {summary_cards}

            <h2>Statistical Summary</h2>
            {stats_table}

            <h2>Prophage Detection</h2>

            <div class="plot-container">
                {prophage_bar}
            </div>

            <div class="plot-container">
                {prophage_box}
            </div>

            <div class="plot-container">
                {detection_tools}
            </div>

            <h2>Quality Assessment</h2>

            <h3>Prophage Quality (CheckV)</h3>
            <div class="plot-container">
                {checkv_quality}
            </div>

            <div class="plot-container">
                {checkv_comp_box}
            </div>

            <h3>MAG Quality (CheckM)</h3>
            <div class="plot-container">
                {checkm_scatter}
            </div>

            <div class="plot-container">
                {checkm_box}
            </div>

            <h2>Taxonomy</h2>

            <div class="plot-container">
                {taxonomy_phylum}
            </div>

            <div class="plot-container">
                {taxonomy_source}
            </div>

            <h2>Assembly Statistics</h2>

            <div class="plot-container">
                {n50_box}
            </div>

            <div class="timestamp">
                Report generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
            </div>
        </div>
    </body>
    </html>
    """

    # Write to file
    with open(output_path, 'w') as f:
        f.write(html)

    print(f"Report generated: {output_path}")

# ============================================================================
# MAIN
# ============================================================================

def main():
    """Main function to generate summary report"""

    # Get inputs from Snakemake
    samples = snakemake.params.samples
    outdir = snakemake.params.outdir
    output_html = snakemake.output.html

    print(f"Generating summary report for {len(samples)} samples...")

    # Parse all data
    print("Parsing prophage tables...")
    prophage_df = parse_prophage_tables(samples, outdir)

    print("Parsing free phage data...")
    free_phage_df = parse_free_phages(samples, outdir)

    print("Parsing CheckV quality data...")
    checkv_all_df = parse_checkv_quality(samples, outdir, 'all_prophages')
    checkv_free_df = parse_checkv_quality(samples, outdir, 'free_phages')

    print("Parsing CheckM quality data...")
    checkm_df = parse_checkm_quality(samples, outdir)

    print("Calculating assembly statistics...")
    assembly_df = parse_assembly_stats(samples, outdir)

    print("Parsing taxonomy data...")
    taxonomy_df = parse_taxonomy(samples, outdir)

    # Generate report
    print("Generating HTML report...")
    generate_html_report(
        prophage_df, free_phage_df, checkv_all_df, checkv_free_df,
        checkm_df, assembly_df, taxonomy_df, output_html
    )

    print("Done!")

if __name__ == '__main__':
    main()
