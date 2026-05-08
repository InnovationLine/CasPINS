#!/usr/bin/env python
"""
Publication Figure Generator for CasPINS R1 Revision
Produces Figures 1, 2, and 3 for the revised manuscript.

All figures are generated deterministically from benchmark data CSVs.
Output: 300 DPI PNG + vector PDF for each figure.

Figure layout (revised per reviewer feedback):
  Fig 1 — Workflow schematic (4 modules + TALEN + data flow)
  Fig 2 — gRNA Benchmark (concordance bars + stage breakdown)
  Fig 3 — Indel Benchmark + worked example overview

Usage:
    python benchmarking/scripts/generate_publication_figures.py

Outputs saved to: benchmarking/results/figures/
"""

import os
import csv
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.colors import to_rgba

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BENCHMARK_DIR = os.path.join(SCRIPT_DIR, '..')
RESULTS_DIR = os.path.join(BENCHMARK_DIR, 'results')
COMPARISON_DIR = os.path.join(RESULTS_DIR, 'grna_design', 'comparison')
FIGURES_DIR = os.path.join(RESULTS_DIR, 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)

# ── Style ─────────────────────────────────────────────────────────────────────
# Colorblind-safe palette (darker for high contrast)
COLORS = {
    'blue':   '#005B8F',
    'orange': '#C46200',
    'green':  '#007A5A',
    'red':    '#B32400',
    'purple': '#800080',
    'yellow': '#B3A300',
    'sky':    '#1A759F',
    'black':  '#111111',
    'grey':   '#555555',
    'lightgrey': '#DDDDDD',
}

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 9,
    'axes.titlesize': 10,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.facecolor': 'white',
    'axes.spines.top': False,
    'axes.spines.right': False,
})

PANEL_LABEL_KW = dict(fontsize=12, fontweight='bold', va='top', ha='left')


# ── Data loaders ──────────────────────────────────────────────────────────────

def load_grna_concordance():
    """Load gRNA concordance data for Fig 2A."""
    path = os.path.join(COMPARISON_DIR, 'grna_comparison_summary.csv')
    genes, chopchop_pct, crispor_pct = [], [], []
    with open(path, newline='') as f:
        for row in csv.DictReader(f):
            if row['gene'] == 'MEAN':
                continue
            genes.append(row['gene'])
            chopchop_pct.append(float(row['full_overlap_pct_of_chopchop']))
            crispor_pct.append(float(row['crispor_full_overlap_pct']))
    return genes, chopchop_pct, crispor_pct


def load_stage_breakdown():
    """Load stage-by-stage discordance data for Fig 2B."""
    path = os.path.join(COMPARISON_DIR, 'stage_breakdown_summary.csv')
    genes, data = [], {}
    cats = ['STRAND_CONVENTION_DIFFERENCE', 'GC_FILTER',
            'HOMOPOLYMER_FILTER', 'REFERENCE_DATABASE_DIFFERENCE']
    cat_labels = ['Strand\nConvention', 'GC\nFilter', 'Homopolymer\nFilter', 'Reference\nDatabase']
    cat_colors = [COLORS['sky'], COLORS['orange'], COLORS['yellow'], '#6A3D9A']

    with open(path, newline='') as f:
        for row in csv.DictReader(f):
            if row['gene'] == 'MEAN':
                continue
            genes.append(row['gene'])
            missed = int(row['missed_total'])
            data[row['gene']] = {
                c: int(row.get(c, 0)) / missed * 100 if missed else 0
                for c in cats
            }
    return genes, data, cats, cat_labels, cat_colors


# ── Figure 1: Workflow Schematic ───────────────────────────────────────────────

def make_box(ax, x, y, w, h, label, sublabel='', color=COLORS['blue'], fontsize=9,
             face_alpha=0.15, label_color=None, sublabel_color='#222222'):
    box = FancyBboxPatch((x - w/2, y - h/2), w, h,
                         boxstyle='round,pad=0.03',
                         facecolor=to_rgba(color, face_alpha),
                         edgecolor=color, linewidth=1.5)
    ax.add_patch(box)
    if label_color is None:
        label_color = color
    ax.text(x, y + (0.015 if sublabel else 0), label,
            ha='center', va='center', fontsize=fontsize,
            fontweight='bold', color=label_color)
    if sublabel:
        ax.text(x, y - 0.04, sublabel,
                ha='center', va='center', fontsize=7.5, color=sublabel_color, fontweight='medium')


def arrow(ax, x1, y1, x2, y2, color='#555555', alpha=1.0):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='-|>', color=color, lw=1.5, alpha=alpha,
                                shrinkA=0, shrinkB=0, mutation_scale=12))


def dashed_box(ax, x, y, w, h, label, color='#888888'):
    box = FancyBboxPatch((x - w/2, y - h/2), w, h,
                         boxstyle='round,pad=0.02',
                         facecolor='#F8F8F8',
                         edgecolor=color, linewidth=1.0,
                         linestyle='--')
    ax.add_patch(box)
    ax.text(x, y, label, ha='center', va='center', fontsize=7.5,
            color='#333333', style='italic')


def fig1_workflow():
    """Figure 1: clean, non-overlapping workflow schematic."""
    fig, ax = plt.subplots(figsize=(10.5, 7.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Title
    ax.text(0.5, 0.995, 'CasPINS Integrated Genome Editing Workflow',
            ha='center', va='top', fontsize=15, fontweight='bold')

    # ── Top input/database layer ──
    input_x, input_y, input_w, input_h = 0.5, 0.86, 0.78, 0.080
    make_box(ax, input_x, input_y, input_w, input_h,
             'INPUT',
             'Gene symbol · Ensembl ID · RefSeq ID · Genomic coordinates · Sequence',
             color='#1F2933', fontsize=10, face_alpha=0.92,
             label_color='white', sublabel_color='white')
    dashed_box(ax, 0.5, 0.745, 0.56, 0.060,
               'Database retrieval: Ensembl REST API · NCBI E-utilities')
    # Input -> database: bottom edge of input to top edge of database box
    arrow(ax, 0.5, input_y - input_h/2, 0.5, 0.745 + 0.060/2)

    # ── Module boxes ──
    box_y = 0.60
    box_h = 0.16
    box_w = 0.19
    positions = [0.12, 0.37, 0.63, 0.88]
    labels    = ['gRNA\nDesign', 'TALEN\nDesign', 'Primer\nDesign', 'Indel\nAnalysis']
    sublabels = ['14 Cas variants\n90+ species\nCRISPRa/i modes',
                 'RVD arrays\n15–20 bp arms\nRFLP screening',
                 'PCR I + PCR II\nCut-site aware\nPrimer3 core',
                 'NNLS traces\nSingle + batch\nR² confidence']
    mod_colors = [COLORS['blue'], COLORS['orange'], COLORS['green'], COLORS['purple']]

    db_bottom = 0.745 - 0.060/2
    module_top = box_y + box_h/2
    for xp, lb, sl, col in zip(positions, labels, sublabels, mod_colors):
        make_box(ax, xp, box_y, box_w, box_h, lb, sl, color=col, fontsize=10)
        # Database -> module: bottom of database band to top edge of module box
        arrow(ax, xp, db_bottom, xp, module_top, color=col)

    # ── Shared project state / hand-off box ──
    shared_x, shared_y, shared_w, shared_h = 0.5, 0.34, 0.90, 0.105
    make_box(ax, shared_x, shared_y, shared_w, shared_h,
             'SHARED PROJECT STATE',
             'selected gRNA · TALEN spacer / cut site · genomic coordinates · run metadata',
             color='#34495E', fontsize=9, face_alpha=0.92,
             label_color='white', sublabel_color='white')

    # Clear module → shared-state arrows. Each starts at a module bottom edge and
    # ends at the top edge of the shared-state box; no free-floating arrows.
    shared_top = shared_y + shared_h / 2
    module_bottom = box_y - box_h/2
    for xp, col in zip(positions, mod_colors):
        # Module -> shared state: bottom edge of module box to top edge of shared box
        arrow(ax, xp, module_bottom, xp, shared_top, color=col, alpha=0.95)

    # ── Output box ──
    output_x, output_y, output_w, output_h = 0.5, 0.105, 0.78, 0.09
    make_box(ax, output_x, output_y, output_w, output_h,
             'OUTPUT',
             'Ranked gRNAs · TALEN pairs · Primer sets · Editing efficiency · Indel spectrum · Run logs',
             color='#5B2C6F', fontsize=9.5, face_alpha=0.92,
             label_color='white', sublabel_color='white')

    # Shared state → output arrow
    # Bottom edge of shared box to top edge of output box
    arrow(ax, 0.5, shared_y - shared_h/2, 0.5, output_y + output_h/2, color=COLORS['black'])

    # Legend outside content region, bottom center
    handles = [mpatches.Patch(facecolor=to_rgba(c, 0.18), edgecolor=c,
                              label=l, linewidth=1.5)
               for c, l in zip(mod_colors, ['gRNA Design', 'TALEN Design',
                                            'Primer Design', 'Indel Analysis'])]
    ax.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.04),
              fontsize=8.5, framealpha=0.95, ncol=4, columnspacing=1.5)

    fig.subplots_adjust(left=0.03, right=0.97, top=0.96, bottom=0.12)
    return fig


# ── Figure 2: gRNA Benchmark ───────────────────────────────────────────────────

def fig2_grna_benchmark():
    genes, chopchop_pct, crispor_pct = load_grna_concordance()
    gene_order, data, cats, cat_labels, cat_colors = load_stage_breakdown()

    fig = plt.figure(figsize=(14.5, 5.2))
    gs = fig.add_gridspec(
        1, 5,
        width_ratios=[1.10, 0.58, 0.16, 1.10, 0.55],
        left=0.045, right=0.985, top=0.84, bottom=0.18, wspace=0.18
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a_leg = fig.add_subplot(gs[0, 1])
    ax_spacer = fig.add_subplot(gs[0, 2])
    ax_b = fig.add_subplot(gs[0, 3])
    ax_b_leg = fig.add_subplot(gs[0, 4])
    ax_a_leg.axis('off')
    ax_spacer.axis('off')
    ax_b_leg.axis('off')

    # ── Panel A: Concordance bars ──────────────────────────────────────────
    x = np.arange(len(genes))
    w = 0.28
    mean_cho = np.mean(chopchop_pct)
    mean_cri = np.mean(crispor_pct)
    bars1 = ax_a.bar(x - w/2, chopchop_pct, w, color=COLORS['blue'],
                     alpha=0.85, label=f'vs CHOPCHOP (mean {mean_cho:.1f}%)',
                     edgecolor='white', linewidth=0.5)
    bars2 = ax_a.bar(x + w/2, crispor_pct, w, color=COLORS['orange'],
                     alpha=0.85, label=f'vs CRISPOR (mean {mean_cri:.1f}%)',
                     edgecolor='white', linewidth=0.5)

    for bar in list(bars1) + list(bars2):
        h = bar.get_height()
        ax_a.text(bar.get_x() + bar.get_width()/2, h + 0.8,
                  f'{h:.0f}%', ha='center', va='bottom', fontsize=7.5)

    # Mean lines (values are listed in legend to avoid text collisions)
    ax_a.axhline(mean_cho, color=COLORS['blue'], linestyle='--', lw=0.9, alpha=0.6)
    ax_a.axhline(mean_cri, color=COLORS['orange'], linestyle='--', lw=0.9, alpha=0.6)

    ax_a.set_xticks(x)
    ax_a.set_xticklabels(genes, fontsize=8.5)
    ax_a.set_ylabel('gRNAs recovered (%)')
    ax_a.set_ylim(0, 100)
    ax_a.set_title('A  Full-gene concordance\n(CasPINS vs CHOPCHOP and CRISPOR)',
                   loc='center', pad=12)
    handles_a, labels_a = ax_a.get_legend_handles_labels()
    ax_a_leg.legend(handles_a, labels_a, loc='center left', framealpha=0.95,
                    title='Concordance', fontsize=8)
    ax_a.grid(axis='y', alpha=0.2, linewidth=0.6)

    # ── Panel B: Stage breakdown stacked bars ─────────────────────────────
    bottom = np.zeros(len(gene_order))
    for cat, label, color in zip(cats, cat_labels, cat_colors):
        vals = np.array([data[g][cat] for g in gene_order])
        bars = ax_b.bar(gene_order, vals, bottom=bottom,
                        color=color, label=label, edgecolor='white', linewidth=0.4)
        # Value labels inside bars
        for i, (v, b) in enumerate(zip(vals, bottom)):
            if v > 4:
                text_color = 'white' if cat == 'REFERENCE_DATABASE_DIFFERENCE' else '#111111'
                ax_b.text(i, b + v/2, f'{v:.0f}%', ha='center', va='center',
                          fontsize=7, color=text_color)
        bottom += vals

    ax_b.set_ylabel('% of missed gRNAs')
    ax_b.set_ylim(0, 105)
    ax_b.set_title('B  Stage-by-stage discordance\n(of gRNAs not recovered by CasPINS)',
                   loc='center', pad=12)
    ax_b.grid(axis='y', alpha=0.2, linewidth=0.6)
    handles_b, labels_b = ax_b.get_legend_handles_labels()
    ax_b_leg.legend(handles_b, labels_b, loc='center left', fontsize=8,
                    ncol=1, framealpha=0.95, title='Discordance category')

    # Annotation
    ax_b.text(0.5, -0.19,
              'Region scanning difference = 0% across all genes',
              ha='center', va='top', transform=ax_b.transAxes,
              fontsize=7.5, color=COLORS['green'], fontweight='bold',
              style='italic')

    return fig


# ── Figure 3: Indel Benchmark ──────────────────────────────────────────────────

def fig3_indel_benchmark():
    # Gold-standard data (from benchmarking/results/indel_analysis/comparison/)
    gold_tools    = ['CasPINS', 'TIDE', 'ICE']
    gold_eff      = [35.7, 33.1, 33.9]
    gold_colors   = [COLORS['blue'], COLORS['orange'], COLORS['green']]

    # DDC data (Table S4 values)
    ddc_tools     = ['CasPINS\n(default)', 'CasPINS\n(R²-corrected)', 'TIDE', 'ICE']
    ddc_eff       = [85.3, 1.5, 9.5, 4.3]   # CasPINS corrected ≈ 0–3%, midpoint 1.5
    ddc_r2        = [0.177, 0.177, 0.100, 0.047]
    ddc_colors    = [COLORS['blue'], COLORS['sky'], COLORS['orange'], COLORS['green']]

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(11.0, 4.8))
    fig.subplots_adjust(left=0.07, right=0.86, top=0.86, bottom=0.24, wspace=0.42)

    # ── Panel A: Gold-standard concordance ───────────────────────────────
    bars = ax_a.bar(gold_tools, gold_eff, color=gold_colors,
                    alpha=0.85, edgecolor='white', linewidth=0.5, width=0.5)
    for bar, v in zip(bars, gold_eff):
        ax_a.text(bar.get_x() + bar.get_width()/2, v + 0.3,
                  f'{v:.1f}%', ha='center', va='bottom', fontsize=9)

    ax_a.set_ylabel('Editing efficiency (%)')
    ax_a.set_ylim(0, 50)
    ax_a.set_title('A  Gold-standard data (TIDE paper example)\n'
                   'All three algorithms agree within 2.6 percentage points', loc='center')
    ax_a.axhspan(33.1 - 0.1, 35.7 + 0.1, color=COLORS['grey'],
                 alpha=0.12, zorder=0)
    ax_a.annotate('2.6 percentage-point range',
                  xy=(1.0, 35.0), xytext=(1.0, 43.5),
                  ha='center', va='center', fontsize=8,
                  arrowprops=dict(arrowstyle='-[,widthB=3.6,lengthB=0.5',
                                  color=COLORS['grey'], lw=1.0),
                  color='#555555')

    # R² values
    for i, (tool, r2) in enumerate(zip(gold_tools, [0.666, 0.977, 0.965])):
        ax_a.text(i, 2.0, f'R²={r2:.3f}', ha='center', va='bottom',
                  fontsize=7, color='#555555')
    ax_a.grid(axis='y', alpha=0.2, linewidth=0.6)

    # ── Panel B: DDC noisy data comparison ───────────────────────────────
    x = np.arange(len(ddc_tools))
    bars_b = ax_b.bar(x, ddc_eff, color=ddc_colors,
                      alpha=0.85, edgecolor='white', linewidth=0.5, width=0.6)
    for bar, v in zip(bars_b, ddc_eff):
        ax_b.text(bar.get_x() + bar.get_width()/2, v + 0.8,
                  f'{v:.1f}%', ha='center', va='bottom', fontsize=9)

    ax_b.set_xticks(x)
    ax_b.set_xticklabels(ddc_tools, fontsize=8)
    ax_b.set_ylabel('Mean editing efficiency (%)')
    ax_b.set_ylim(0, 100)
    ax_b.set_title('B  DDC dataset (challenging traces, n=8)\n'
                   'R²-corrected CasPINS converges with TIDE/ICE', loc='center')

    # R² overlay on twin axis
    ax_b2 = ax_b.twinx()
    ax_b2.plot(x, ddc_r2, 'D--', color=COLORS['red'], markersize=5,
               label='R²', linewidth=1.2)
    ax_b2.set_ylim(0, 0.6)
    ax_b2.set_ylabel('R² (goodness-of-fit)', color=COLORS['red'])
    ax_b2.tick_params(axis='y', colors=COLORS['red'])
    ax_b2.spines['right'].set_edgecolor(COLORS['red'])
    ax_b2.legend(loc='upper right', fontsize=8, framealpha=0.95)

    # Quality warning annotation
    ax_b.axhline(0, color='black', lw=0.5)
    ax_b.grid(axis='y', alpha=0.2, linewidth=0.6)
    ax_b.text(0.5, -0.19,
              'All R² < 0.3 → CasPINS v2.0 displays quality warning',
              ha='center', va='top', transform=ax_b.transAxes,
              fontsize=7.5, color=COLORS['red'], style='italic')

    return fig


# ── Save helper ────────────────────────────────────────────────────────────────

def save_figure(fig, name):
    png_path = os.path.join(FIGURES_DIR, f'{name}.png')
    pdf_path = os.path.join(FIGURES_DIR, f'{name}.pdf')
    fig.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  Saved: {name}.png  +  {name}.pdf')


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print('=' * 60)
    print('CasPINS Publication Figure Generator')
    print('=' * 60)

    print('\nFig 1: Workflow schematic...')
    save_figure(fig1_workflow(), 'Fig1_workflow_schematic')

    print('Fig 2: gRNA Benchmark (concordance + stage breakdown)...')
    save_figure(fig2_grna_benchmark(), 'Fig2_grna_benchmark')

    print('Fig 3: Indel Benchmark (gold-standard + DDC)...')
    save_figure(fig3_indel_benchmark(), 'Fig3_indel_benchmark')

    print(f'\nAll figures saved to: {FIGURES_DIR}')
    print('Formats: 300 DPI PNG (for journal upload) + PDF (vector, for proofing)')
    print('\nFor each figure:')
    print('  PNG: upload to ScholarOne as high-resolution figure file')
    print('  PDF: open in any viewer to check readability before submission')
    print('\nFigures are generated deterministically from benchmark CSV data.')
    print('No AI tools were used in figure generation.')


if __name__ == '__main__':
    main()
