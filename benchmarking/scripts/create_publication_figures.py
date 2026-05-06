#!/usr/bin/env python
"""
Create publication-quality figures for Bioinformatics Advances submission.

This script:
1. Reformats manuscript figures to 600 DPI TIFF with standardized panel labels
2. Converts supplementary tables (S1-S4) to publication figure format
3. Follows OUP guidelines: 600 DPI, TIFF with LZW compression, centered

Usage:
    python benchmarking/scripts/create_publication_figures.py
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont
import numpy as np
import csv
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BENCHMARK_DIR = os.path.join(SCRIPT_DIR, '..')
OUT_DIR = os.path.join(BENCHMARK_DIR, 'results', 'manuscript_figures', 'publication')
EXTRACTED_DIR = os.path.join(BENCHMARK_DIR, 'results', 'manuscript_figures', 'extracted')
SUPP_DIR = os.path.join(BENCHMARK_DIR, 'results', 'supplementary')
os.makedirs(OUT_DIR, exist_ok=True)

DPI = 600
FULL_PAGE_WIDTH_INCHES = 7.0
TARGET_WIDTH_PX = int(FULL_PAGE_WIDTH_INCHES * DPI)
PANEL_FONT_SIZE = 72

try:
    font_bold = ImageFont.truetype("arialbd.ttf", PANEL_FONT_SIZE)
except Exception:
    try:
        font_bold = ImageFont.truetype("arial.ttf", PANEL_FONT_SIZE)
    except Exception:
        font_bold = ImageFont.load_default()


# ================================================================
# PART 1: Manuscript Figure Reformatting
# ================================================================

FIGURE_MAP = {
    'Fig1_gRNA_Design': {
        'src': os.path.join(EXTRACTED_DIR, 'figure_6.jpeg'),
        'panels': ['A', 'B'],
    },
    'Fig2_TALEN_Design': {
        'src': os.path.join(EXTRACTED_DIR, 'figure_2.jpeg'),
        'panels': ['A', 'B', 'C'],
    },
    'Fig3_Primer_Design_Input': {
        'src': os.path.join(EXTRACTED_DIR, 'figure_4.jpeg'),
        'panels': ['A', 'B'],
    },
    'Fig3_Primer_Design_Results': {
        'src': os.path.join(EXTRACTED_DIR, 'figure_7.jpeg'),
        'panels': ['C', 'D'],
    },
    'Fig4_Single_Indel': {
        'src': os.path.join(EXTRACTED_DIR, 'figure_5.jpeg'),
        'panels': ['A', 'B', 'C'],
    },
    'Fig5_Batch_Analysis': {
        'src': os.path.join(EXTRACTED_DIR, 'figure_1.jpeg'),
        'panels': ['A', 'B', 'C'],
    },
    'Graphical_Abstract': {
        'src': os.path.join(EXTRACTED_DIR, 'figure_3.jpg'),
        'panels': [],
    },
}


def process_manuscript_figure(name, info):
    src = info['src']
    if not os.path.exists(src):
        print(f"  SKIP: {src} not found")
        return

    img = Image.open(src).convert('RGB')
    orig_w, orig_h = img.size
    scale = TARGET_WIDTH_PX / orig_w
    new_w = TARGET_WIDTH_PX
    new_h = int(orig_h * scale)
    img = img.resize((new_w, new_h), Image.LANCZOS)

    draw = ImageDraw.Draw(img)
    panels = info['panels']
    margin = 20

    if panels:
        n = len(panels)
        positions = [(margin, i * new_h // n + margin) for i in range(n)]

        for letter, (x, y) in zip(panels, positions):
            try:
                bbox = font_bold.getbbox(letter)
                tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
            except Exception:
                tw, th = PANEL_FONT_SIZE, PANEL_FONT_SIZE

            pad = 10
            draw.rectangle([x - pad, y - pad, x + tw + pad * 2, y + th + pad * 2],
                           fill='white')
            draw.text((x, y), letter, fill='black', font=font_bold)

    out_path = os.path.join(OUT_DIR, f'{name}.tiff')
    img.save(out_path, format='TIFF', dpi=(DPI, DPI), compression='tiff_lzw')
    print(f"  {name}: {new_w}x{new_h} px, {DPI} DPI -> {out_path}")


def process_all_manuscript_figures():
    print("=" * 60)
    print("PART 1: MANUSCRIPT FIGURE REFORMATTING")
    print(f"Target: {DPI} DPI, {FULL_PAGE_WIDTH_INCHES}\" width, TIFF/LZW")
    print("=" * 60)
    for name, info in FIGURE_MAP.items():
        process_manuscript_figure(name, info)
    print()


# ================================================================
# PART 2: Supplementary Table Figures
# ================================================================

def create_table_s1():
    """Table S1: Workflow Timing Comparison"""
    print("  Creating Table S1: Workflow Timing Comparison...")

    data = [
        ['gRNA Design (per gene)', 'CRISPOR', '5-8', '1-2',
         'Requires NCBI sequence fetch + paste'],
        ['gRNA Design (per gene)', 'CHOPCHOP', '7-10', '1-2',
         'Longer computation times'],
        ['Indel Analysis (per sample)', 'TIDE', '2-3', '0.5-1',
         'Separate control + edited uploads'],
        ['Indel Analysis (per sample)', 'ICE', '3-5', '0.5-1',
         'Requires template spreadsheet'],
        ['Primer Design (per target)', 'Primer-BLAST', '5-8', '0.5-1',
         'Manual sequence extraction'],
        ['Data Collation', 'Manual', '5-10', '0',
         'CasPINS auto-integrates results'],
        ['Total (5 genes + 8 samples)', 'Multiple', '79-143', '9-19',
         '3-5 separate web tools required'],
    ]

    col_labels = ['Task', 'Traditional Tool', 'Traditional\nTime (min)',
                  'CasPINS\nTime (min)', 'Notes']

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.axis('off')
    ax.set_title('Supplementary Table S1: Workflow Timing Comparison\n'
                 'CasPINS Integrated Workflow vs. Traditional Fragmented Workflow',
                 fontsize=11, fontweight='bold', pad=20, loc='left')

    table = ax.table(cellText=data, colLabels=col_labels, loc='center',
                     cellLoc='center', colColours=['#2c5282'] * 5)
    table.auto_set_font_size(False)
    table.set_fontsize(8)

    for key, cell in table.get_celld().items():
        row, col = key
        if row == 0:
            cell.set_text_props(color='white', fontweight='bold', fontsize=8.5)
            cell.set_facecolor('#2c5282')
            cell.set_height(0.08)
        else:
            if row == len(data):
                cell.set_text_props(fontweight='bold')
                cell.set_facecolor('#e8f0fe')
            elif row % 2 == 0:
                cell.set_facecolor('#f7fafc')
            else:
                cell.set_facecolor('white')
            cell.set_height(0.08)
        cell.set_edgecolor('#cbd5e0')
        cell.set_linewidth(0.5)

    widths = [0.22, 0.13, 0.12, 0.12, 0.30]
    for ci, w in enumerate(widths):
        for ri in range(len(data) + 1):
            table[ri, ci].set_width(w)

    fig.text(0.05, 0.04,
             'Time estimates based on standard user workflow with modern hardware and broadband internet.\n'
             'CasPINS times assume the GUI is pre-loaded. Bold row = total for a typical experiment.',
             fontsize=6.5, style='italic', color='#4a5568')

    out = os.path.join(OUT_DIR, 'Table_S1_workflow_timing.tiff')
    fig.savefig(out, dpi=DPI, bbox_inches='tight', format='tiff',
                pil_kwargs={'compression': 'tiff_lzw'})
    plt.close()
    print(f"    -> {out}")


def create_table_s2():
    """Table S2: gRNA Comparison Summary"""
    print("  Creating Table S2: gRNA Comparison Summary...")

    data = [
        ['ATE1', 'ENSG00000107669', '37,911', '129', '190', '72.1', '74.5', '60.5', '0.029'],
        ['DBH', 'ENSG00000123454', '6,189', '457', '334', '64.1', '67.1', '59.5', '-0.015'],
        ['EMX1', 'ENSG00000170370', '5,019', '206', '289', '65.0', '65.7', '49.5', '-0.094'],
        ['TP53', 'ENSG00000141510', '5,781', '105', '307', '74.3', '66.8', '45.7', '-0.120'],
        ['VEGFA', 'ENSG00000112715', '4,151', '60', '356', '68.3', '61.8', '48.3', '-0.129'],
        ['Mean', '', '11,810', '191', '295', '68.8', '67.2', '52.7', '-0.066'],
    ]

    col_labels = ['Gene', 'Ensembl ID', 'CasPINS\ngRNAs', 'CHOPCHOP\ngRNAs', 'CRISPOR\ngRNAs',
                  'vs CHOPCHOP\nOverlap (%)', 'vs CRISPOR\nOverlap (%)',
                  'Region\nOverlap (%)', 'Spearman\n\u03C1']

    fig, ax = plt.subplots(figsize=(11, 4.5))
    ax.axis('off')
    ax.set_title('Supplementary Table S2: gRNA Design Comparison Summary\n'
                 'CasPINS vs. CHOPCHOP and CRISPOR across 5 Human Genes (hg38, SpCas9 NGG)',
                 fontsize=10, fontweight='bold', pad=20, loc='left')

    table = ax.table(cellText=data, colLabels=col_labels, loc='center',
                     cellLoc='center', colColours=['#2c5282'] * 9)
    table.auto_set_font_size(False)
    table.set_fontsize(8)

    for key, cell in table.get_celld().items():
        row, col = key
        if row == 0:
            cell.set_text_props(color='white', fontweight='bold', fontsize=7.5)
            cell.set_facecolor('#2c5282')
            cell.set_height(0.12)
        else:
            if row == len(data):
                cell.set_text_props(fontweight='bold')
                cell.set_facecolor('#e8f0fe')
            elif row % 2 == 0:
                cell.set_facecolor('#f7fafc')
            else:
                cell.set_facecolor('white')
            cell.set_height(0.09)
        cell.set_edgecolor('#cbd5e0')
        cell.set_linewidth(0.5)

    widths = [0.06, 0.15, 0.08, 0.08, 0.08, 0.11, 0.11, 0.09, 0.08]
    for ci, w in enumerate(widths):
        for ri in range(len(data) + 1):
            table[ri, ci].set_width(w)

    fig.text(0.05, 0.06,
             'Full-gene overlap: percentage of external tool gRNAs found in CasPINS full-gene output.\n'
             'Region overlap: CasPINS filtered to CHOPCHOP target region. '
             'Spearman \u03C1 computed on shared gRNAs.',
             fontsize=6.5, style='italic', color='#4a5568')

    out = os.path.join(OUT_DIR, 'Table_S2_grna_comparison.tiff')
    fig.savefig(out, dpi=DPI, bbox_inches='tight', format='tiff',
                pil_kwargs={'compression': 'tiff_lzw'})
    plt.close()
    print(f"    -> {out}")


def create_table_s3():
    """Table S3: Top 20 gRNAs per Gene"""
    print("  Creating Table S3: Top 20 gRNAs per Gene...")

    rows = []
    csv_path = os.path.join(SUPP_DIR, 'Table_S3_top20_grnas_per_gene.csv')
    with open(csv_path) as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            if row:
                rows.append(row)

    genes = ['TP53', 'ATE1', 'VEGFA', 'DBH', 'EMX1']

    fig, axes = plt.subplots(5, 1, figsize=(11, 16))
    fig.suptitle('Supplementary Table S3: Top 10 gRNAs per Gene from CasPINS\n'
                 '(Human, hg38, SpCas9 NGG)',
                 fontsize=11, fontweight='bold', y=0.99)

    gene_colors = {
        'TP53': '#2c5282', 'ATE1': '#276749',
        'VEGFA': '#9b2c2c', 'DBH': '#553c9a', 'EMX1': '#744210'
    }

    for gene_idx, gene in enumerate(genes):
        ax = axes[gene_idx]
        ax.axis('off')

        gene_rows = [r for r in rows if r[0] == gene][:10]

        table_data = []
        for r in gene_rows:
            table_data.append([
                r[1], r[2], r[3], r[4],
                f"{float(r[5]):.0f}%", f"{float(r[7]):.4f}",
                f"{float(r[8]):.3f}", f"{float(r[9]):.3f}", f"{float(r[10]):.3f}",
            ])

        col_labels = ['Rank', 'gRNA Sequence (20nt)', 'PAM', 'Strand',
                      'GC%', 'Composite', 'Doench', 'Moreno-M.', 'Xu']

        color = gene_colors.get(gene, '#2c5282')
        ax.set_title(f'{gene}', fontsize=9, fontweight='bold', loc='left', pad=3)

        table = ax.table(cellText=table_data, colLabels=col_labels, loc='center',
                         cellLoc='center', colColours=[color] * 9)
        table.auto_set_font_size(False)
        table.set_fontsize(5.5)

        for key, cell in table.get_celld().items():
            r, c = key
            if r == 0:
                cell.set_text_props(color='white', fontweight='bold', fontsize=5.5)
                cell.set_facecolor(color)
                cell.set_height(0.035)
            else:
                cell.set_facecolor('#f7fafc' if r % 2 == 0 else 'white')
                cell.set_height(0.028)
            cell.set_edgecolor('#e2e8f0')
            cell.set_linewidth(0.3)

        widths = [0.04, 0.25, 0.05, 0.06, 0.05, 0.08, 0.07, 0.07, 0.07]
        for ci, w in enumerate(widths):
            for ri in range(len(table_data) + 1):
                table[ri, ci].set_width(w)

    fig.text(0.05, 0.005,
             'Top 10 of 20 gRNAs shown per gene. Composite Score = weighted average of '
             'Doench 2016, Moreno-Mateos, and Xu on-target efficiency scores. '
             'Full 20 per gene in CSV.',
             fontsize=5.5, style='italic', color='#4a5568')

    plt.subplots_adjust(hspace=0.45, top=0.96)

    out = os.path.join(OUT_DIR, 'Table_S3_top20_grnas.tiff')
    fig.savefig(out, dpi=DPI, bbox_inches='tight', format='tiff',
                pil_kwargs={'compression': 'tiff_lzw'})
    plt.close()
    print(f"    -> {out}")


def create_table_s4():
    """Table S4: Indel Analysis Comparison"""
    print("  Creating Table S4: Indel Analysis Comparison...")

    data = [
        ['Open Source', 'GFP-target', 'example2', '35.7', '0.666', '33.1', '0.977', '33.9', '0.965'],
        ['DDC', 'DDC', 'editedA2', '100.0', '0.276', '14.5', '0.157', '12.9', '0.145'],
        ['DDC', 'DDC', 'editedA3', '70.7', '0.263', '10.4', '0.113', '4.4', '0.045'],
        ['DDC', 'DDC', 'editedA4', '67.8', '0.130', '3.1', '0.038', '0.2', '0.002'],
        ['DDC', 'DDC', 'editedA7', '94.5', '0.289', '14.4', '0.146', '5.0', '0.054'],
        ['DDC', 'DDC', 'editedB6', '89.8', '0.398', '18.2', '0.182', '4.5', '0.046'],
        ['DDC', 'DDC', 'editedB7', '87.6', '0.000', '2.0', '0.021', '1.4', '0.015'],
        ['DDC', 'DDC', 'editedB10', '85.3', '0.000', '3.8', '0.041', '1.8', '0.021'],
        ['DDC', 'DDC', 'editedB12', '86.8', '0.063', '9.9', '0.105', '4.4', '0.048'],
        ['DDC', 'DDC', 'MEAN', '85.3', '0.177', '9.5', '0.100', '4.3', '0.047'],
    ]

    col_labels = ['Dataset', 'Gene', 'Sample',
                  'CasPINS\nEff. (%)', 'CasPINS\nR\u00B2',
                  'TIDE\nEff. (%)', 'TIDE\nR\u00B2',
                  'ICE\nEff. (%)', 'ICE\nR\u00B2']

    fig, ax = plt.subplots(figsize=(11, 6.5))
    ax.axis('off')
    ax.set_title('Supplementary Table S4: Three-Way Indel Analysis Algorithm Comparison\n'
                 'CasPINS vs. TIDE (Brinkman et al. 2014) vs. ICE (Hsiau et al. 2019)',
                 fontsize=10, fontweight='bold', pad=20, loc='left')

    table = ax.table(cellText=data, colLabels=col_labels, loc='center',
                     cellLoc='center', colColours=['#553c9a'] * 9)
    table.auto_set_font_size(False)
    table.set_fontsize(7.5)

    for key, cell in table.get_celld().items():
        row, col = key
        if row == 0:
            cell.set_text_props(color='white', fontweight='bold', fontsize=7.5)
            cell.set_facecolor('#553c9a')
            cell.set_height(0.07)
        elif row == 1:
            cell.set_facecolor('#e9d8fd')
            cell.set_text_props(fontweight='bold')
            cell.set_height(0.05)
        elif row == len(data):
            cell.set_text_props(fontweight='bold')
            cell.set_facecolor('#feebc8')
            cell.set_height(0.05)
        else:
            cell.set_facecolor('#faf5ff' if row % 2 == 0 else 'white')
            cell.set_height(0.05)
        cell.set_edgecolor('#d6bcfa')
        cell.set_linewidth(0.5)

    widths = [0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09]
    for ci, w in enumerate(widths):
        for ri in range(len(data) + 1):
            table[ri, ci].set_width(w)

    fig.text(0.05, 0.06,
             'Algorithms implemented faithfully from published methods:\n'
             'TIDE: Stacked 4-channel peak heights, single NNLS, R\u00B2-corrected '
             '(Brinkman et al. 2014, NAR 42:e168)\n'
             'ICE: Normalized peaks, Lasso (L1, \u03B1=0.8), R\u00B2-corrected '
             '(Hsiau et al. 2019, CRISPR J 2:123-130)\n'
             'CasPINS: NNLS on combined (summed) trace channels',
             fontsize=5.5, style='italic', color='#4a5568')
    fig.text(0.05, 0.02,
             'Purple row = Open source gold-standard data: TIDE (Brinkman et al. 2014, gku936 Supplementary) vs ICE (Hsiau et al. 2019). '
             'Orange row = DDC experimental mean. '
             'DDC control has only 5 basecalled bases.',
             fontsize=5.5, style='italic', color='#4a5568')

    out = os.path.join(OUT_DIR, 'Table_S4_indel_comparison.tiff')
    fig.savefig(out, dpi=DPI, bbox_inches='tight', format='tiff',
                pil_kwargs={'compression': 'tiff_lzw'})
    plt.close()
    print(f"    -> {out}")


def create_all_table_figures():
    print("\n" + "=" * 60)
    print("PART 2: SUPPLEMENTARY TABLE FIGURES")
    print(f"Target: {DPI} DPI, TIFF with LZW compression")
    print("=" * 60)
    create_table_s1()
    create_table_s2()
    create_table_s3()
    create_table_s4()
    print()


if __name__ == '__main__':
    process_all_manuscript_figures()
    create_all_table_figures()
    print("=" * 60)
    print("ALL FIGURES COMPLETE")
    print(f"Output directory: {OUT_DIR}")
    print("=" * 60)
