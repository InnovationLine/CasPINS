#!/usr/bin/env python
"""
CasPINS - Cas-Primer-Indel Suite - GUI Interface

A user-friendly graphical interface for CRISPR gRNA design, primer design, and indel analysis.
This GUI wraps the existing command-line tools while preserving all core functionality.

Enhanced Features:
- 90+ species support (CHOPCHOP-level)
- Multiple editing modes: Knockout, Knock-in (HDR), CRISPRa, CRISPRi
- TALEN support
- RefSeq/Ensembl/Genomic coordinates input
- TIDE-style visualization
- In-GUI results display
"""

import streamlit as st
import os
import sys
import subprocess
import json
import pandas as pd
from datetime import datetime
import tempfile
import shutil
from pathlib import Path
import logging
import threading

# Try to import tkinter for folder picker
try:
    import tkinter as tk
    from tkinter import filedialog
    TKINTER_AVAILABLE = True
except ImportError:
    TKINTER_AVAILABLE = False

# Configure logging for debugging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)


def convert_species_to_ensembl(species_input: str) -> str:
    """Convert species display name to Ensembl format (genus_species)."""
    # Comprehensive mapping from display names/common names to Ensembl species names
    species_map = {
        # Common names to Ensembl format
        'human': 'homo_sapiens',
        'homo sapiens': 'homo_sapiens',
        'mouse': 'mus_musculus',
        'mus musculus': 'mus_musculus',
        'rat': 'rattus_norvegicus',
        'rattus norvegicus': 'rattus_norvegicus',
        'cat': 'felis_catus',
        'felis catus': 'felis_catus',
        'dog': 'canis_lupus_familiaris',
        'canis familiaris': 'canis_lupus_familiaris',
        'canis lupus familiaris': 'canis_lupus_familiaris',
        'zebrafish': 'danio_rerio',
        'danio rerio': 'danio_rerio',
        'chicken': 'gallus_gallus',
        'gallus gallus': 'gallus_gallus',
        'pig': 'sus_scrofa',
        'sus scrofa': 'sus_scrofa',
        'cow': 'bos_taurus',
        'cattle': 'bos_taurus',
        'bos taurus': 'bos_taurus',
        'sheep': 'ovis_aries',
        'ovis aries': 'ovis_aries',
        'goat': 'capra_hircus',
        'capra hircus': 'capra_hircus',
        'horse': 'equus_caballus',
        'equus caballus': 'equus_caballus',
        'rabbit': 'oryctolagus_cuniculus',
        'oryctolagus cuniculus': 'oryctolagus_cuniculus',
        'guinea pig': 'cavia_porcellus',
        'cavia porcellus': 'cavia_porcellus',
        'hamster': 'mesocricetus_auratus',
        'mesocricetus auratus': 'mesocricetus_auratus',
        'frog': 'xenopus_tropicalis',
        'xenopus': 'xenopus_tropicalis',
        'xenopus tropicalis': 'xenopus_tropicalis',
        'fruit fly': 'drosophila_melanogaster',
        'drosophila': 'drosophila_melanogaster',
        'drosophila melanogaster': 'drosophila_melanogaster',
        'worm': 'caenorhabditis_elegans',
        'c. elegans': 'caenorhabditis_elegans',
        'caenorhabditis elegans': 'caenorhabditis_elegans',
        'yeast': 'saccharomyces_cerevisiae',
        'saccharomyces cerevisiae': 'saccharomyces_cerevisiae',
        'macaque': 'macaca_mulatta',
        'rhesus': 'macaca_mulatta',
        'macaca mulatta': 'macaca_mulatta',
        'cynomolgus': 'macaca_fascicularis',
        'macaca fascicularis': 'macaca_fascicularis',
        'chimpanzee': 'pan_troglodytes',
        'pan troglodytes': 'pan_troglodytes',
        'gorilla': 'gorilla_gorilla',
        'gorilla gorilla': 'gorilla_gorilla',
        'orangutan': 'pongo_abelii',
        'pongo abelii': 'pongo_abelii',
        'marmoset': 'callithrix_jacchus',
        'callithrix jacchus': 'callithrix_jacchus',
        'bonobo': 'pan_paniscus',
        'pan paniscus': 'pan_paniscus',
        'medaka': 'oryzias_latipes',
        'oryzias latipes': 'oryzias_latipes',
        'pufferfish': 'takifugu_rubripes',
        'fugu': 'takifugu_rubripes',
        'takifugu rubripes': 'takifugu_rubripes',
        'stickleback': 'gasterosteus_aculeatus',
        'gasterosteus aculeatus': 'gasterosteus_aculeatus',
        'salmon': 'salmo_salar',
        'salmo salar': 'salmo_salar',
        'tilapia': 'oreochromis_niloticus',
        'oreochromis niloticus': 'oreochromis_niloticus',
        'turkey': 'meleagris_gallopavo',
        'meleagris gallopavo': 'meleagris_gallopavo',
        'duck': 'anas_platyrhynchos',
        'anas platyrhynchos': 'anas_platyrhynchos',
        'platypus': 'ornithorhynchus_anatinus',
        'ornithorhynchus anatinus': 'ornithorhynchus_anatinus',
        'opossum': 'monodelphis_domestica',
        'monodelphis domestica': 'monodelphis_domestica',
        'elephant': 'loxodonta_africana',
        'loxodonta africana': 'loxodonta_africana',
        'ferret': 'mustela_putorius_furo',
        'mustela putorius furo': 'mustela_putorius_furo',
        'armadillo': 'dasypus_novemcinctus',
        'dasypus novemcinctus': 'dasypus_novemcinctus',
        'alpaca': 'vicugna_pacos',
        'vicugna pacos': 'vicugna_pacos',
        'dolphin': 'tursiops_truncatus',
        'tursiops truncatus': 'tursiops_truncatus',
    }
    
    # Clean up input - extract from display format "Human (Homo sapiens)"
    species_clean = species_input.split('(')[0].strip().lower()
    
    # Also check parenthetical content
    if '(' in species_input and ')' in species_input:
        paren_content = species_input.split('(')[1].split(')')[0].strip().lower()
    else:
        paren_content = ''
    
    # Try exact match first
    if species_clean in species_map:
        return species_map[species_clean]
    
    if paren_content and paren_content in species_map:
        return species_map[paren_content]
    
    # Try partial match
    for key, value in species_map.items():
        if key in species_clean or species_clean in key:
            return value
        if paren_content and (key in paren_content or paren_content in key):
            return value
    
    # If no match, try to convert to Ensembl format (genus_species)
    # Remove extra words and convert spaces to underscores
    cleaned = species_clean.replace(' ', '_').replace('-', '_')
    
    # Log warning for unknown species
    logger.warning(f"Unknown species mapping: '{species_input}' -> trying '{cleaned}'")
    
    return cleaned

# Add parent (src) directory to path for imports
_src_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _src_dir not in sys.path:
    sys.path.insert(0, _src_dir)

# Import existing modules
from grna_design import GRNADesigner
from grna_design.grna_designer import EditingMode, TargetType
from grna_design.hdr_designer import HDRDesigner, EditType
from grna_design.talen_designer import TALENDesigner

from utils.primer_design import (
    generate_primer_recommendations, 
    generate_enhanced_primer_report,
    PrimerDesigner
)
from utils.sequence_analysis import find_grna_in_sequence
from utils.run_logger import RunLogger
from utils.ab1_parser import parse_ab1
from utils.indel_analysis import decompose_traces_indel_analysis, calculate_editing_efficiency_fallback
from utils.visualization_multi import plot_indel_analysis_multi, create_summary_report
from utils.file_management import ensure_output_dir, archive_output_files
from utils.tide_visualization import (
    create_tide_style_figure,
    display_results_streamlit,
    display_indel_spectrum_streamlit,
    display_indel_spectrum_simple,
    figure_to_base64
)

# Import species configuration
try:
    from config import get_species_for_dropdown, get_crispr_systems, get_editing_modes, get_species_by_name
except ImportError:
    # Fallback if config not available
    def get_species_for_dropdown():
        return [
            "Human (Homo sapiens)", "Mouse (Mus musculus)", "Rat (Rattus norvegicus)",
            "Zebrafish (Danio rerio)", "Fruit fly (Drosophila melanogaster)",
            "Roundworm (Caenorhabditis elegans)", "Chicken (Gallus gallus)",
            "Pig (Sus scrofa)", "Cow (Bos taurus)", "Dog (Canis lupus familiaris)"
        ]
    def get_crispr_systems():
        return [{"name": "SpCas9", "pam": "NGG"}, {"name": "SaCas9", "pam": "NNGRRT"},
                {"name": "Cas12a", "pam": "TTTV"}]
    def get_editing_modes():
        return [{"name": "knockout", "description": "Gene disruption via NHEJ"}]
    def get_species_by_name(name):
        return None

# Import settings configuration
try:
    from config.settings import (
        get_data_directory, set_data_directory, is_data_directory_configured,
        list_available_genes, get_gene_folder, load_config, save_config
    )
except ImportError:
    # Fallback - will use session state for data directory
    def get_data_directory():
        return st.session_state.get('data_directory', '')
    def set_data_directory(path):
        st.session_state['data_directory'] = path
        return True
    def is_data_directory_configured():
        return bool(st.session_state.get('data_directory', ''))
    def list_available_genes():
        return []
    def get_gene_folder(gene_name):
        data_dir = get_data_directory()
        if data_dir:
            folder = os.path.join(data_dir, gene_name.lower())
            return folder if os.path.isdir(folder) else None
        return None
    def load_config():
        return {}
    def save_config(config):
        pass

# Set page configuration
st.set_page_config(
    page_title="CasPINS - Cas-Primer-Indel Suite",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
   .main {
        padding-top: 1rem;
    }
   .stTabs [data-baseweb="tab-list"] {
        gap: 24px;
    }
   .stTabs [data-baseweb="tab"] {
        height: 50px;
        padding-left: 20px;
        padding-right: 20px;
    }
   .success-box {
        padding: 1rem;
        border-radius: 0.5rem;
        background-color: #d4edda;
        border: 1px solid #c3e6cb;
        color: #155724;
        margin: 1rem 0;
    }
   .info-box {
        padding: 1rem;
        border-radius: 0.5rem;
        background-color: #d1ecf1;
        border: 1px solid #bee5eb;
        color: #0c5460;
        margin: 1rem 0;
    }
   .warning-box {
        padding: 1rem;
        border-radius: 0.5rem;
        background-color: #fff3cd;
        border: 1px solid #ffeaa7;
        color: #856404;
        margin: 1rem 0;
    }
</style>
""", unsafe_allow_html=True)

# Folder picker function using tkinter
def open_folder_picker():
    """Open a native folder picker dialog and return the selected path."""
    if not TKINTER_AVAILABLE:
        return None
    
    selected_path = None
    
    def run_dialog():
        nonlocal selected_path
        root = tk.Tk()
        root.withdraw()  # Hide the main window
        root.attributes('-topmost', True)  # Bring dialog to front
        
        # Open folder selection dialog
        selected_path = filedialog.askdirectory(
            title="Select Data Directory for CRISPR Analysis",
            initialdir=os.path.expanduser("~")
        )
        
        root.destroy()
    
    # Run in main thread (required for tkinter)
    run_dialog()
    
    return selected_path if selected_path else None


# Initialize session state
if 'results' not in st.session_state:
    st.session_state.results = {}
if 'grna_results' not in st.session_state:
    st.session_state.grna_results = None
if 'selected_grnas' not in st.session_state:
    st.session_state.selected_grnas = []
if 'grna_full_data' not in st.session_state:
    st.session_state.grna_full_data = {}
if 'folder_picker_result' not in st.session_state:
    st.session_state.folder_picker_result = None

def main():
    """Main GUI application."""

    # Initialize session state for data directory
    if 'data_directory' not in st.session_state:
        st.session_state.data_directory = get_data_directory() or ""

    # Header
    st.title("🧬 CasPINS - Cas-Primer-Indel Suite")

    # Sidebar for general settings
    with st.sidebar:
        st.header("⚙️ Settings")

        # Data directory management (CRITICAL)
        st.subheader("📁 Data Directory")
        st.markdown("*Select folder for your gene data*")
        
        # Current data directory display
        current_data_dir = st.session_state.data_directory
        
        if current_data_dir and os.path.isdir(current_data_dir):
            st.success(f"📂 **Active:** `{current_data_dir}`")
        else:
            st.warning("⚠️ No data directory set")
        
        # BROWSE BUTTON - Primary way to select folder
        st.markdown("---")
        
        if TKINTER_AVAILABLE:
            if st.button("📂 Browse & Select Folder", type="primary", width="stretch", key="browse_folder_btn"):
                selected_folder = open_folder_picker()
                if selected_folder:
                    st.session_state.data_directory = selected_folder
                    set_data_directory(selected_folder)
                    st.rerun()
            st.caption("Click to open folder picker dialog")
        else:
            st.info("📝 Enter path manually below")
        
        # Manual path entry (alternative)
        with st.expander("✏️ Or Enter Path Manually", expanded=not TKINTER_AVAILABLE):
            new_data_dir = st.text_input(
                "Data Folder Path",
                value=current_data_dir,
                placeholder="C:/Users/YourName/crispr_data",
                help="Full path to directory containing gene folders",
                key="sidebar_data_dir"
            )
            
            col1, col2 = st.columns(2)
            with col1:
                if st.button("📂 Set Path", key="set_data_dir"):
                    if new_data_dir:
                        if os.path.isdir(new_data_dir):
                            st.session_state.data_directory = new_data_dir
                            set_data_directory(new_data_dir)
                            st.success("✅ Set!")
                            st.rerun()
                        else:
                            st.error("❌ Not found!")
                    else:
                        st.warning("⚠️ Enter a path")
            
            with col2:
                if st.button("📝 Create", key="create_data_dir"):
                    if new_data_dir:
                        try:
                            os.makedirs(new_data_dir, exist_ok=True)
                            st.session_state.data_directory = new_data_dir
                            set_data_directory(new_data_dir)
                            st.success("✅ Created!")
                            st.rerun()
                        except Exception as e:
                            st.error(f"❌ {e}")
        
        # Quick set options
        with st.expander("⚡ Quick Options", expanded=False):
            quick_paths = [
                ("Home/crispr_data", os.path.join(os.path.expanduser("~"), "crispr_data")),
                ("Documents/crispr_data", os.path.join(os.path.expanduser("~"), "Documents", "crispr_data")),
                ("Current/data", os.path.join(os.getcwd(), "data")),
            ]
            
            for label, qpath in quick_paths:
                if st.button(f"📁 {label}", key=f"quick_{hash(qpath)}", width="stretch"):
                    try:
                        os.makedirs(qpath, exist_ok=True)
                        st.session_state.data_directory = qpath
                        set_data_directory(qpath)
                        st.rerun()
                    except Exception as e:
                        st.error(f"❌ {e}")
        
        # Show status and Project Structure
        data_dir = st.session_state.data_directory
        if data_dir and os.path.isdir(data_dir):
            st.markdown("---")
            st.subheader("📊 Project Structure")
            
            with st.expander("ℹ️ How data is stored", expanded=False):
                st.markdown("""
                CasPINS organizes data by gene:
                ```text
                data_directory/
                └── gene_name/
                    ├── grna.txt
                    ├── mrna.txt
                    ├── primers/
                    ├── input/  (for indels)
                    │   ├── control.ab1
                    │   └── edited.ab1
                    └── output/ (results)
                ```
                """)
            
            # List existing genes and their contents
            existing_genes = [d for d in os.listdir(data_dir)
                            if os.path.isdir(os.path.join(data_dir, d)) and not d.startswith('.')]
            
            if existing_genes:
                with st.expander("📁 Current Project Files", expanded=True):
                    for gene in sorted(existing_genes):
                        gene_path = os.path.join(data_dir, gene)
                        
                        # Count files inside
                        grna_exists = "✅" if os.path.exists(os.path.join(gene_path, "grna.txt")) else "❌"
                        mrna_exists = "✅" if os.path.exists(os.path.join(gene_path, "mrna.txt")) else "❌"
                        
                        # Check inputs
                        input_path = os.path.join(gene_path, "input")
                        ab1_count = len([f for f in os.listdir(input_path) if f.endswith('.ab1')]) if os.path.exists(input_path) else 0
                        
                        st.markdown(f"**{gene.upper()}**")
                        st.caption(f"gRNA: {grna_exists} | mRNA: {mrna_exists} | AB1 files: {ab1_count}")
            
            # Create new gene folder
            with st.expander("➕ Create New Gene Folder"):
                new_gene = st.text_input("Gene Name", placeholder="e.g., MYOD1", key="sidebar_new_gene")
                if st.button("Create", key="sidebar_create_folder") and new_gene:
                    gene_folder = os.path.join(data_dir, new_gene.lower())
                    os.makedirs(gene_folder, exist_ok=True)
                    os.makedirs(os.path.join(gene_folder, "input"), exist_ok=True)
                    # Create template files
                    with open(os.path.join(gene_folder, "grna.txt"), 'w') as f:
                        f.write("# Add gRNA sequences here (one per line)\n")
                    with open(os.path.join(gene_folder, "mrna.txt"), 'w') as f:
                        f.write("# Add mRNA/cDNA sequence here\n")
                    st.success(f"✅ Created: {gene_folder}")
                    st.rerun()
        else:
            st.warning("⚠️ Please set a valid data directory")
            st.info("💡 Required for Indel Analysis and local file gRNA design")
    
    # Main tabs
    tab1, tab2, tab3, tab4 = st.tabs(["🔍 gRNA Design", "🧪 Primer Design", "📊 Indel Analysis", "📚 Documentation"])
    
    with tab1:
        grna_design_tab()
    
    with tab2:
        primer_design_tab()
    
    with tab3:
        indel_analysis_tab()
    
    with tab4:
        documentation_tab()


def grna_design_tab():
    """gRNA Design interface with enhanced features."""
    st.header("🔍 gRNA Design")
    st.markdown("Find optimal guide RNAs for your target gene - supports **90+ species**, multiple input types, and editing modes")
    
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("Target Selection")
        
        # Input type selection
        input_type = st.radio(
            "Input Type",
            ["Gene Symbol", "Ensembl ID", "RefSeq ID", "Genomic Coordinates", "Paste Sequence"],
            horizontal=True,
            key="grna_input_type"
        )
        
        # Gene input based on type
        if input_type == "Gene Symbol":
            gene_name = st.text_input("Gene Symbol", 
                                     placeholder="e.g., TP53, DDC, BRCA1",
                                     help="Enter official gene symbol",
                                     key="grna_gene_input")
        elif input_type == "Ensembl ID":
            gene_name = st.text_input("Ensembl Gene/Transcript ID", 
                                     placeholder="e.g., ENSG00000107669 or ENST00000358428",
                                     help="Enter Ensembl gene or transcript ID",
                                     key="grna_ensembl_input")
        elif input_type == "RefSeq ID":
            gene_name = st.text_input("RefSeq Accession", 
                                     placeholder="e.g., NM_003054.5 or NG_008849.1",
                                     help="Enter RefSeq mRNA or gene accession",
                                     key="grna_refseq_input")
        elif input_type == "Genomic Coordinates":
            gene_name = st.text_input("Genomic Coordinates", 
                                     placeholder="e.g., chr8:20100000-20150000",
                                     help="Enter coordinates as chr:start-end",
                                     key="grna_coord_input")
        else:  # Paste Sequence
            gene_name = st.text_area("DNA Sequence", 
                                    placeholder="Paste your DNA sequence here (FASTA or raw)",
                                    height=100,
                                    key="grna_seq_input")
        
        # Species selection - Expanded list (90+ species)
        species_options = get_species_for_dropdown()
        species_display = st.selectbox(
            "Species", 
            species_options,
            index=0, 
            key="grna_species_expanded",
            help="Select from 90+ supported species"
        )
        # Extract simple species name
        species = species_display.split('(')[0].strip().lower().replace(' ', '_')
        
        # Get assemblies for selected species
        species_info = get_species_by_name(species_display.split('(')[0].strip())
        assemblies = []
        if species_info and 'assemblies' in species_info:
            assemblies = [f"{a['name']} ({a.get('alias', '')})" if a.get('alias') else a['name'] for a in species_info['assemblies']]
            
        assembly_options = ["Auto (Default)"] + assemblies

        # Assembly Selection
        assembly_selection = st.selectbox(
            "Genome Assembly",
            assembly_options,
            index=0,
            help="Select specific assembly for this species, or use the default",
            key="grna_assembly"
        )
        
        if assembly_selection == "Auto (Default)":
            assembly = None
        else:
            assembly = assembly_selection.split(' ')[0].strip()
        
        # Editing Mode - NEW
        st.markdown("#### Editing Mode")
        editing_mode = st.selectbox(
            "Application",
            [
                "Knockout (NHEJ)",
                "Knock-in (HDR)",
                "CRISPRa (Activation)",
                "CRISPRi (Repression)",
                "Base Editing (C>T)",
                "Base Editing (A>G)",
                "Nanopore Enrichment"
            ],
            index=0,
            key="grna_editing_mode",
            help="Select the type of editing experiment"
        )
        
        # Show editing mode info
        mode_info = {
            "Knockout (NHEJ)": "Target early exons for gene disruption",
            "Knock-in (HDR)": "Design gRNAs near insertion site (<10bp from desired edit)",
            "CRISPRa (Activation)": "Target -200 to -50 bp upstream of TSS",
            "CRISPRi (Repression)": "Target -50 to +300 bp around TSS",
            "Base Editing (C>T)": "C must be in positions 4-8 of protospacer",
            "Base Editing (A>G)": "A must be in positions 4-7 of protospacer",
            "Nanopore Enrichment": "Design flanking gRNAs for region enrichment"
        }
        st.info(f"💡 {mode_info.get(editing_mode, '')}")
        
        # CRISPR system vs TALEN
        nuclease_type = st.radio("Nuclease System", ["CRISPR-Cas", "TALEN"], horizontal=True, key="grna_nuclease_type")
        
        if nuclease_type == "CRISPR-Cas":
            cas_type = st.selectbox("CRISPR-Cas System",
                                   ["SpCas9 (NGG)", "SpCas9-NG (NG)", "SaCas9 (NNGRRT)", 
                                    "Cas12a/Cpf1 (TTTV)", "xCas9 (NG/GAA/GAT)",
                                    "SpCas9-VQR (NGAN)", "SpCas9-EQR (NGAG)", 
                                    "SpCas9-VRER (NGCG)", "CjCas9 (NNNNRYAC)"],
                                   key="grna_cas_type_expanded")
            # Extract cas type name
            cas_type_simple = cas_type.split('(')[0].strip().replace('-', '_').replace('/', '_')
        else:
            st.info("🔧 TALEN design uses paired binding sites with RVD sequences")
            cas_type_simple = "TALEN"
    
    with col2:
        st.subheader("Filters & Scoring")
        
        # Filters
        st.markdown("#### Sequence Filters")
        with st.expander("Advanced Filters", expanded=True):
            col2a, col2b = st.columns(2)
            with col2a:
                gc_min = st.slider("Min GC%", 0, 100, 40, key="grna_gc_min")
                homopolymer_max = st.number_input("Max Homopolymer", min_value=3, max_value=8, value=4, key="grna_homopolymer_max")
            with col2b:
                gc_max = st.slider("Max GC%", 0, 100, 60, key="grna_gc_max")
                top_n = st.number_input("Number of gRNAs", min_value=1, max_value=50, value=10, key="grna_top_n")
            
            no_poly_t = st.checkbox("Remove poly-T in seed region", key="grna_no_poly_t")
            check_off_targets = st.checkbox("Check off-targets (slower)", key="grna_check_off_targets")
        
        # Target region (only for knockout mode)
        if "Knockout" in editing_mode:
            target_region = st.selectbox("Target Region",
                                        ["Early exons (recommended)", "All exons", "CDS only", 
                                         "5' UTR", "3' UTR", "All regions"],
                                        key="grna_target_region_new")
        else:
            target_region = "all"
        
        # Display current settings
        st.markdown("#### Current Settings Summary")
        settings_data = {
            'Parameter': ['Species', 'Nuclease', 'Mode', 'GC Range', 'Results'],
            'Value': [species_display.split('(')[0].strip(), 
                     cas_type_simple if nuclease_type == "CRISPR-Cas" else "TALEN",
                     editing_mode.split('(')[0].strip(), 
                     f"{gc_min}-{gc_max}%", 
                     str(top_n)]
        }
        settings_df = pd.DataFrame(settings_data)
        st.dataframe(settings_df, hide_index=True, width="stretch")
    
    # Design buttons
    col_btn1, col_btn2 = st.columns([3, 1])
    
    with col_btn1:
        find_clicked = st.button("🚀 Find gRNAs", type="primary", width="stretch", key="grna_find_button")
    
    with col_btn2:
        if st.button("🗑️ Clear Results", width="stretch", key="grna_clear_button"):
            st.session_state.grna_results = None
            st.session_state.selected_grnas = []
            st.session_state.grna_full_data = {}
            st.rerun()
    
    if find_clicked:
        if not gene_name:
            st.error("Please enter a gene name or Ensembl ID")
            return
        
        # Branch based on nuclease type
        if nuclease_type == "TALEN":
            # TALEN Design
            with st.spinner(f"Searching for TALEN pairs in {gene_name}..."):
                try:
                    logger.info(f"TALEN Design: target={gene_name}, species={species}")
                    logger.info(f"Input type: {input_type}")
                    
                    # Initialize TALEN designer
                    talen_designer = TALENDesigner(
                        species=species,
                        assembly=assembly if assembly else None
                    )
                    
                    # Convert species to Ensembl format for display
                    ensembl_species = convert_species_to_ensembl(species_display)
                    
                    # Show search information
                    st.info(f"🔍 Searching for TALEN sites in: **{gene_name}** ({ensembl_species})")
                    if assembly:
                        st.info(f"📊 Assembly: {assembly}")
                    
                    # Determine target type based on input
                    if input_type == "Gene Symbol":
                        target_type = 'gene'
                    elif input_type == "Ensembl ID":
                        target_type = 'ensembl_gene'
                    else:
                        target_type = 'sequence'
                        
                    logger.info(f"Target type resolved to: {target_type}")
                    
                    data_dir = st.session_state.get('data_directory', '')
                    run_log = RunLogger(data_dir, "talen_design", gene_name)
                    run_log.log_parameters({
                        "species": species, "assembly": assembly,
                        "n_results": top_n, "target_type": target_type,
                        "check_off_targets": check_off_targets
                    })

                    # Design TALENs
                    results = talen_designer.design_talens(
                        target=gene_name,
                        n_results=top_n,
                        target_type=target_type,
                        include_off_targets=check_off_targets
                    )
                    
                    num_pairs = len(results.get('talen_pairs', []))
                    logger.info(f"TALEN Design results: {num_pairs} pairs found")
                    run_log.log_results({"talen_pairs_found": num_pairs})
                    
                    if 'error' in results:
                        logger.warning(f"Design error: {results['error']}")
                        run_log.log_error(results['error'])
                        
                    log_file_path = run_log.finish()
                    st.success(f"Log saved: {log_file_path}")
                    
                    # Store results (mark as TALEN results)
                    results['nuclease_type'] = 'TALEN'
                    st.session_state.grna_results = results
                    
                    # Check for errors
                    if 'error' in results:
                        st.error(f"Error: {results['error']}")
                        st.info("💡 **Troubleshooting:** Try using an Ensembl ID instead of gene symbol")
                        return
                    
                    # Results will be displayed outside the button block
                    talen_pairs = results.get('talen_pairs', [])
                    if talen_pairs:
                        st.success(f"Found {len(talen_pairs)} TALEN pairs! (See results below)")
                    else:
                        st.warning("No suitable TALEN pairs found")
                        st.info("The sequence may not have valid TALEN target sites")
                        
                except Exception as e:
                    st.error(f"Error during TALEN design: {str(e)}")
                    logger.exception("TALEN design error")
                    st.info("Please check your input and try again")
        else:
            # CRISPR-Cas Design
            with st.spinner(f"Searching for {cas_type_simple} gRNAs in {gene_name}..."):
                try:
                    # Get data directory for local file lookup
                    data_dir = st.session_state.get('data_directory', '')
                    
                    logger.info(f"gRNA Design: target={gene_name}, species={species}, cas={cas_type_simple}")
                    logger.info(f"Input type: {input_type}, Data dir: {data_dir or 'Not set'}")
                    
                    # Initialize designer - fetches from online databases (like CHOPCHOP)
                    designer = GRNADesigner(
                        species=species,
                        assembly=assembly if assembly else None,
                        cas_type=cas_type_simple
                    )
                    
                    # Prepare filters
                    filters = {
                        'gc_min': gc_min,
                        'gc_max': gc_max,
                        'homopolymer_max': homopolymer_max,
                        'remove_poly_t_seed': no_poly_t
                    }
                    
                    # Convert species to Ensembl format for display
                    ensembl_species = convert_species_to_ensembl(species_display)
                    
                    # Show search information
                    st.info(f"🔍 Searching for: **{gene_name}** in **{ensembl_species}**")
                    if assembly:
                        st.info(f"📊 Assembly: {assembly}")
                    
                    # Determine target type based on input
                    if input_type == "Gene Symbol":
                        target_type = 'gene'
                    elif input_type == "Ensembl ID":
                        target_type = 'ensembl_gene'
                        # Warn if Ensembl ID might not match selected species
                        if gene_name.upper().startswith('ENSG') and 'homo_sapiens' not in ensembl_species:
                            st.warning(f"⚠️ Note: ENSG IDs are human-specific. Looking up gene symbol to search in {ensembl_species}.")
                        elif gene_name.upper().startswith('ENSRNOG') and 'rattus' not in ensembl_species:
                            st.warning(f"⚠️ Note: ENSRNOG IDs are rat-specific. Looking up gene symbol to search in {ensembl_species}.")
                        elif gene_name.upper().startswith('ENSMUSG') and 'mus_musculus' not in ensembl_species:
                            st.warning(f"⚠️ Note: ENSMUSG IDs are mouse-specific. Looking up gene symbol to search in {ensembl_species}.")
                    elif input_type == "RefSeq ID":
                        target_type = 'refseq_mrna'
                    elif input_type == "Genomic Coordinates":
                        target_type = 'genomic_coords'
                    else:
                        target_type = 'sequence'
                    
                    logger.info(f"Target type resolved to: {target_type}")
                    
                    data_dir = st.session_state.get('data_directory', '')
                    run_log = RunLogger(data_dir, "grna_design", gene_name)
                    run_log.log_parameters({
                        "species": species, "assembly": assembly, "cas_type": cas_type_simple,
                        "n_results": top_n, "filters": filters, "target_type": target_type,
                        "check_off_targets": check_off_targets
                    })
                    
                    # Design gRNAs
                    results = designer.design_grnas(
                        target=gene_name,
                        n_results=top_n,
                        filters=filters,
                        include_off_targets=check_off_targets,
                        target_type=target_type
                    )
                    
                    num_grnas = len(results.get('grnas', []))
                    logger.info(f"Design results: {num_grnas} gRNAs found")
                    run_log.log_results({"grnas_found": num_grnas})
                    
                    if 'error' in results:
                        logger.warning(f"Design error: {results['error']}")
                        run_log.log_error(results['error'])
                        
                    log_file_path = run_log.finish()
                    st.success(f"Log saved: {log_file_path}")
                    
                    # Store results
                    results['nuclease_type'] = 'CRISPR'
                    st.session_state.grna_results = results
                    
                    # Check for errors
                    if 'error' in results:
                        st.error(f"Error: {results['error']}")
                        
                        # Provide specific troubleshooting based on the error
                        if "Could not retrieve sequence" in results['error']:
                            st.info("💡 **Troubleshooting:**")
                            st.markdown("""
                            1. **Check gene name spelling** - Try the official gene symbol
                            2. **Try Ensembl ID** - More reliable than gene names (e.g., ENSRNOG00000011218)
                            3. **Check species** - Make sure the species is correct
                            4. **Use local file** - Create `<data_dir>/{}/mrna.txt` with your sequence
                            """.format(gene_name.lower()))
                            
                            # Check if local file exists and suggest creating it
                            data_dir = st.session_state.get('data_directory', '')
                            local_file = os.path.join(data_dir, gene_name.lower(), "mrna.txt") if data_dir else ""
                            if not os.path.exists(local_file):
                                st.markdown(f"**Create local file:** `{local_file}` with your mRNA sequence")
                        
                        return
                    
                    # Results will be displayed outside the button block
                    grnas = results.get('grnas', [])
                    if grnas:
                        st.success(f"Found {len(grnas)} high-quality gRNAs! (See results below)")
                    else:
                        st.warning("No suitable gRNAs found with current filters")
                        st.info("Try adjusting GC content range or other filters")
                        
                except Exception as e:
                    st.error(f"Error during gRNA design: {str(e)}")
                    st.info("Please check your input and try again")
    
    # Display results OUTSIDE the button block so they persist on checkbox clicks
    if st.session_state.grna_results:
        results = st.session_state.grna_results
        nuclease_type_result = results.get('nuclease_type', 'CRISPR')
        
        if nuclease_type_result == 'TALEN' and 'talen_pairs' in results:
            # Display TALEN results
            talen_pairs = results.get('talen_pairs', [])
            if talen_pairs:
                st.markdown("---")
                st.header("🔧 TALEN Design Results")
                
                # Display gene metadata
                metadata = results.get('metadata', {})
                if metadata:
                    with st.expander("📋 Gene Information", expanded=True):
                        col1, col2 = st.columns(2)
                        with col1:
                            st.markdown(f"**Gene:** {metadata.get('gene_name', 'N/A')}")
                            st.markdown(f"**Ensembl ID:** {metadata.get('ensembl_id', 'N/A')}")
                            st.markdown(f"**Species:** {metadata.get('species', 'N/A')}")
                        with col2:
                            st.markdown(f"**Chromosome:** {metadata.get('chromosome', 'N/A')}")
                            st.markdown(f"**Assembly:** {metadata.get('assembly', 'N/A')}")
                            st.markdown(f"**Sequence Length:** {metadata.get('sequence_length', 'N/A')} bp")
                
                # Display TALEN pair results
                display_talen_results(talen_pairs, metadata)
        
        elif 'grnas' in results:
            # Display CRISPR gRNA results
            grnas = results.get('grnas', [])
            
            if grnas:
                st.markdown("---")
                st.header("🧬 gRNA Design Results")
                
                # Display gene metadata
                metadata = results.get('metadata', {})
                if metadata:
                    with st.expander("📋 Gene Information", expanded=True):
                        col1, col2 = st.columns(2)
                        with col1:
                            st.markdown(f"**Gene:** {metadata.get('gene_name', 'N/A')}")
                            st.markdown(f"**Symbol:** {metadata.get('gene_symbol', 'N/A')}")
                            st.markdown(f"**Ensembl ID:** {metadata.get('ensembl_id', 'N/A')}")
                        with col2:
                            st.markdown(f"**Chromosome:** {metadata.get('chromosome', 'N/A')}")
                            strand_display = '+ (sense)' if metadata.get('strand') == 1 else '- (antisense)'
                            st.markdown(f"**Strand:** {strand_display}")
                            st.markdown(f"**Assembly:** {metadata.get('assembly', 'N/A')}")
                
                # Display gRNAs in a table
                display_grna_results(grnas)


def display_talen_results(talen_pairs: list, metadata: dict = None):
    """
    Display TALEN design results with color-coded sequences and comprehensive analysis.
    
    Features:
    - Color-coded target sequence (TALE1 blue, spacer gray, TALE2 green)
    - Multi-row RVD display format
    - Off-target analysis
    - Restriction site positions
    - Separate Rank (availability) vs Best ID (efficiency)
    """
    st.subheader("TALEN Pair Results")
    
    # Information about TALEN design
    with st.expander("ℹ️ About TALEN Design", expanded=False):
        st.markdown("""
        **TALEN (Transcription Activator-Like Effector Nuclease)** pairs consist of:
        - **Left arm (TALE 1)**: Binds 5' of the spacer on the sense strand
        - **Right arm (TALE 2)**: Binds 5' of the spacer on the antisense strand  
        - **Spacer region**: Where the FokI nuclease domains create a double-strand break
        
        **RVD Codes**: NI=A, HD=C, NG=T, NN/NH=G
        
        **Optimal parameters**:
        - Arm length: 15-20 bp (optimal: 17-18 bp)
        - Spacer length: 12-21 bp (optimal: 14-16 bp)
        - GC content: 40-60%
        
        **Ranking**:
        - **Rank**: Sequential order based on target position/availability
        - **Best ID**: Ranking based on predicted editing efficiency
        """)
    
    def format_colored_sequence(pair):
        """Create color-coded HTML for target sequence."""
        left_seq = pair.get('left_arm_seq', '')
        spacer_seq = pair.get('spacer_seq', '')
        right_seq = pair.get('right_arm_seq', '')
        
        # Colors: TALE1 = blue, Spacer = gray, TALE2 = green
        html = f'<span style="color: #1E88E5; font-family: monospace; font-weight: bold;">{left_seq}</span>'
        html += f'<span style="color: #757575; font-family: monospace;">{spacer_seq}</span>'
        html += f'<span style="color: #43A047; font-family: monospace; font-weight: bold;">{right_seq}</span>'
        return html
    
    def format_rvd_multiline(rvd_string, cols=6):
        """Format RVD sequence into multi-row HTML table."""
        if not rvd_string:
            return "-"
        rvds = rvd_string.replace(' ', '-').split('-')
        rows = []
        for i in range(0, len(rvds), cols):
            row = rvds[i:i+cols]
            rows.append(' '.join(row))
        return '<br>'.join(rows)
    
    def format_restriction_sites(sites):
        """Format restriction sites list for display."""
        if not sites or sites == '-':
            return '-'
        if isinstance(sites, list):
            return ' '.join(sites)
        return str(sites)
    
    # Summary table
    st.markdown("### Results Summary")
    
    # Prepare data for table
    table_data = []
    for pair in talen_pairs:
        # Format restriction sites
        sites = pair.get('restriction_sites', [])
        sites_display = format_restriction_sites(sites)
        
        # Format target sequence (truncate if too long for display)
        target_seq = pair.get('target_sequence', '')
        if len(target_seq) > 55:
            target_display = target_seq[:25] + '...' + target_seq[-25:]
        else:
            target_display = target_seq
        
        row = {
            'Rank': pair.get('rank', 0),
            'Target Sequence': target_display,
            'Genomic Location': pair.get('genomic_location', 'N/A'),
            'Cluster': pair.get('cluster', 0),
            'OT Pairs': pair.get('off_target_pairs', 0),
            # Use formatted display strings (found/total) like CHOPCHOP
            'MM0': pair.get('mm0_display', f"{pair.get('off_targets_mm0', 0)}/{pair.get('off_targets_mm0_total', 1)}"),
            'MM1': pair.get('mm1_display', f"{pair.get('off_targets_mm1', 0)}/{pair.get('off_targets_mm1_total', 10)}"),
            'MM2': pair.get('mm2_display', f"{pair.get('off_targets_mm2', 0)}/{pair.get('off_targets_mm2_total', 50)}"),
            'MM3': pair.get('mm3_display', f"{pair.get('off_targets_mm3', 0)}/{pair.get('off_targets_mm3_total', 80)}"),
            'Restriction Sites': sites_display,
            'Best ID': pair.get('best_id', 0),
            'Efficiency': f"{pair.get('efficiency_score', 0):.2f}",
        }
        table_data.append(row)
    
    # Create DataFrame
    df = pd.DataFrame(table_data)
    
    # Display the table
    st.dataframe(
        df,
        width="stretch",
        hide_index=True,
        column_config={
            "Rank": st.column_config.NumberColumn("Rank", help="Sequential rank by position", width="small"),
            "Target Sequence": st.column_config.TextColumn("Target Sequence", help="Full target sequence (TALE1+Spacer+TALE2)", width="large"),
            "Genomic Location": st.column_config.TextColumn("Genomic Location", width="medium"),
            "Cluster": st.column_config.NumberColumn("Cluster", help="Nearby TALEN sites", width="small"),
            "OT Pairs": st.column_config.NumberColumn("OT Pairs", help="Off-target pairs", width="small"),
            "MM0": st.column_config.TextColumn("MM0", help="0 mismatches (found/total)", width="small"),
            "MM1": st.column_config.TextColumn("MM1", help="1 mismatch (found/total)", width="small"),
            "MM2": st.column_config.TextColumn("MM2", help="2 mismatches (found/total)", width="small"),
            "MM3": st.column_config.TextColumn("MM3", help="3 mismatches (found/total)", width="small"),
            "Restriction Sites": st.column_config.TextColumn("Restriction Sites", help="Enzyme:Position", width="medium"),
            "Best ID": st.column_config.NumberColumn("Best ID", help="Rank by efficiency", width="small"),
            "Efficiency": st.column_config.TextColumn("Efficiency", help="Predicted editing efficiency (0-1)", width="small"),
        }
    )
    
    # Detailed view for each TALEN pair
    st.markdown("### Detailed TALEN Pair Information")
    
    for i, pair in enumerate(talen_pairs):
        rank = pair.get('rank', i+1)
        best_id = pair.get('best_id', i+1)
        efficiency = pair.get('efficiency_score', 0)
        
        with st.expander(f"TALEN Pair #{rank} (Best ID: {best_id}) - Efficiency: {efficiency:.3f}", expanded=(i == 0)):
            # Color-coded target sequence
            st.markdown("**Target Sequence** (color-coded: "
                       "<span style='color: #1E88E5;'>TALE1</span> | "
                       "<span style='color: #757575;'>Spacer</span> | "
                       "<span style='color: #43A047;'>TALE2</span>)", 
                       unsafe_allow_html=True)
            
            colored_seq = format_colored_sequence(pair)
            st.markdown(f"<div style='background: #f5f5f5; padding: 10px; border-radius: 5px; overflow-x: auto;'>{colored_seq}</div>", 
                       unsafe_allow_html=True)
            
            st.markdown("---")
            
            # TALE1 and TALE2 RVD sequences in multi-row format
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("**TALE 1 (Left Arm)**")
                st.markdown(f"**Sequence:** `{pair.get('left_arm_seq', 'N/A')}`")
                st.markdown(f"**Length:** {pair.get('left_length', len(pair.get('left_arm_seq', '')))} bp")
                st.markdown(f"**GC Content:** {pair.get('left_gc', 0):.1f}%")
                
                # Multi-row RVD display
                rvd_html = format_rvd_multiline(pair.get('tale1_raw', pair.get('tale1', '')))
                st.markdown("**RVD Sequence:**")
                st.markdown(f"<div style='background: #E3F2FD; padding: 8px; border-radius: 4px; font-family: monospace;'>{rvd_html}</div>", 
                           unsafe_allow_html=True)
            
            with col2:
                st.markdown("**TALE 2 (Right Arm)**")
                st.markdown(f"**Sequence:** `{pair.get('right_arm_seq', 'N/A')}`")
                st.markdown(f"**Length:** {pair.get('right_length', len(pair.get('right_arm_seq', '')))} bp")
                st.markdown(f"**GC Content:** {pair.get('right_gc', 0):.1f}%")
                
                # Multi-row RVD display
                rvd_html = format_rvd_multiline(pair.get('tale2_raw', pair.get('tale2', '')))
                st.markdown("**RVD Sequence:**")
                st.markdown(f"<div style='background: #E8F5E9; padding: 8px; border-radius: 4px; font-family: monospace;'>{rvd_html}</div>", 
                           unsafe_allow_html=True)
            
            # Spacer region
            st.markdown("**Spacer Region**")
            spacer_col1, spacer_col2 = st.columns([3, 1])
            with spacer_col1:
                st.code(pair.get('spacer_seq', 'N/A'), language=None)
            with spacer_col2:
                st.metric("Length", f"{pair.get('spacer_length', 0)} bp")
            
            # Off-target Analysis (displayed as found/total like CHOPCHOP)
            st.markdown("**Off-target Analysis** (found/total sites)")
            ot_cols = st.columns(6)
            with ot_cols[0]:
                st.metric("OT Pairs", pair.get('off_target_pairs', 0))
            with ot_cols[1]:
                mm0_display = pair.get('mm0_display', f"{pair.get('off_targets_mm0', 0)}/{pair.get('off_targets_mm0_total', 1)}")
                st.metric("MM0", mm0_display, help="0 mismatches (exact)")
            with ot_cols[2]:
                mm1_display = pair.get('mm1_display', f"{pair.get('off_targets_mm1', 0)}/{pair.get('off_targets_mm1_total', 10)}")
                st.metric("MM1", mm1_display, help="1 mismatch")
            with ot_cols[3]:
                mm2_display = pair.get('mm2_display', f"{pair.get('off_targets_mm2', 0)}/{pair.get('off_targets_mm2_total', 50)}")
                st.metric("MM2", mm2_display, help="2 mismatches")
            with ot_cols[4]:
                mm3_display = pair.get('mm3_display', f"{pair.get('off_targets_mm3', 0)}/{pair.get('off_targets_mm3_total', 80)}")
                st.metric("MM3", mm3_display, help="3 mismatches")
            with ot_cols[5]:
                specificity = pair.get('specificity_score', 0)
                st.metric("Specificity", f"{specificity:.2f}", help="Higher is better")
            
            # Restriction sites
            sites = pair.get('restriction_sites', [])
            if sites and sites != '-':
                st.markdown("**Restriction Sites** (Enzyme:Position)")
                sites_display = format_restriction_sites(sites)
                st.info(sites_display)
            else:
                st.markdown("**Restriction Sites:** None found")
            
            # Genomic location and scores
            score_cols = st.columns(3)
            with score_cols[0]:
                if pair.get('genomic_location') and pair.get('genomic_location') != 'N/A':
                    st.markdown(f"**Genomic Location:** `{pair.get('genomic_location')}`")
            with score_cols[1]:
                st.markdown(f"**Composite Score:** {pair.get('composite_score', 0):.3f}")
            with score_cols[2]:
                st.markdown(f"**Cluster Size:** {pair.get('cluster', 0)}")
    
    # Download options
    st.markdown("### 📥 Download Results")
    
    col1, col2 = st.columns(2)
    
    with col1:
        # CSV download
        csv_data = df.to_csv(index=False)
        st.download_button(
            label="Download Summary (CSV)",
            data=csv_data,
            file_name=f"talen_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
            mime="text/csv"
        )
    
    with col2:
        # Detailed JSON download
        json_data = json.dumps(talen_pairs, indent=2, default=str)
        st.download_button(
            label="Download Full Details (JSON)",
            data=json_data,
            file_name=f"talen_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
            mime="application/json"
        )


def display_grna_results(grnas, metadata=None):
    """Display gRNA results in an interactive table with genomic locations."""
    st.subheader("gRNA Results")
    
    # Show genomic context if available
    if metadata:
        col1, col2, col3 = st.columns(3)
        with col1:
            if metadata.get('genomic_location'):
                st.metric("Genomic Location", metadata['genomic_location'])
        with col2:
            if metadata.get('ensembl_id'):
                st.metric("Ensembl ID", metadata['ensembl_id'])
        with col3:
            if metadata.get('refseq_id'):
                st.metric("RefSeq ID", metadata['refseq_id'])
    
    # Store full gRNA data for detail view
    if 'grna_full_data' not in st.session_state:
        st.session_state.grna_full_data = {}
    
    # Prepare data for display with enhanced information
    data = []
    for i, grna in enumerate(grnas):
        # Store full data keyed by sequence
        st.session_state.grna_full_data[str(grna['sequence'])] = grna
        
        row = {
            'Select': False,
            'Rank': i + 1,
            'Sequence': str(grna['sequence']),
            'PAM': str(grna.get('pam', 'N/A')),
            'Position': str(grna.get('position', 'N/A')),
            'Strand': str(grna.get('strand', 'N/A')),
            'GC%': f"{grna.get('gc_content', 0):.1f}",
            'Score': f"{grna['scores']['composite']:.3f}",
        }
        
        # Add genomic location if available (CHOPCHOP-style)
        if grna.get('genomic_location'):
            row['Genomic Location'] = grna['genomic_location']
        elif grna.get('location_simple'):
            row['Genomic Location'] = grna['location_simple']
        else:
            row['Genomic Location'] = 'N/A'
        
        # Add detailed scores
        row['Doench'] = f"{grna['scores'].get('doench_2016', 0):.3f}"
        row['MM'] = f"{grna['scores'].get('moreno_mateos', 0):.3f}"
        row['Xu'] = f"{grna['scores'].get('xu', 0):.3f}"
        
        # Add base editing info if present
        if grna.get('base_editing_targets'):
            row['Edit Positions'] = str(grna['base_editing_targets'])
        
        data.append(row)
    
    # Create DataFrame
    df = pd.DataFrame(data)
    
    # Configure columns based on what's available
    column_config = {
        "Select": st.column_config.CheckboxColumn(
            "Select",
            help="Select gRNAs to view details and save",
            default=False,
        ),
        "Sequence": st.column_config.TextColumn(
            "gRNA Sequence",
            help="20nt guide sequence",
            width="medium",
        ),
        "Genomic Location": st.column_config.TextColumn(
            "Genomic Location",
            help="Chromosomal position (chr:start-end)",
            width="medium",
        ),
        "Score": st.column_config.NumberColumn(
            "Overall Score",
            help="Composite score (0-1)",
            format="%.3f",
        ),
    }
    
    # Columns to disable
    disabled_cols = ["Rank", "Sequence", "PAM", "Position", "Strand", "GC%", 
                    "Score", "Doench", "MM", "Xu", "Genomic Location"]
    if 'Edit Positions' in df.columns:
        disabled_cols.append('Edit Positions')
    
    # Display with selection
    edited_df = st.data_editor(
        df,
        hide_index=True,
        width="stretch",
        column_config=column_config,
        disabled=disabled_cols,
        key="grna_selector"
    )
    
    # Get selected gRNAs
    selected_rows = edited_df[edited_df['Select']]
    selected_grnas = selected_rows['Sequence'].tolist()
    st.session_state.selected_grnas = selected_grnas
    
    if selected_grnas:
        st.success(f"✅ Selected {len(selected_grnas)} gRNA(s)")
        
        # Show detailed view for selected gRNAs
        st.markdown("---")
        st.subheader("📋 Selected gRNA Details")
        
        for seq in selected_grnas:
            grna_data = st.session_state.grna_full_data.get(seq, {})
            
            with st.expander(f"🧬 gRNA: {seq}", expanded=True):
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    st.markdown("**Basic Information**")
                    st.markdown(f"- **Sequence:** `{seq}`")
                    st.markdown(f"- **PAM:** `{grna_data.get('pam', 'NGG')}`")
                    st.markdown(f"- **Full Target:** `{seq}{grna_data.get('pam', 'NGG')}`")
                    st.markdown(f"- **Length:** {len(seq)} nt")
                    st.markdown(f"- **GC Content:** {grna_data.get('gc_content', 0):.1f}%")
                    st.markdown(f"- **Strand:** {grna_data.get('strand', 'N/A')}")
                
                with col2:
                    st.markdown("**Genomic Position**")
                    st.markdown(f"- **Position in Gene:** {grna_data.get('position', 'N/A')}")
                    if grna_data.get('genomic_location'):
                        st.markdown(f"- **Genomic Location:** `{grna_data['genomic_location']}`")
                    elif grna_data.get('location_simple'):
                        st.markdown(f"- **Location:** `{grna_data['location_simple']}`")
                    if grna_data.get('exon'):
                        st.markdown(f"- **Exon:** {grna_data['exon']}")
                
                with col3:
                    st.markdown("**Scoring Details**")
                    scores = grna_data.get('scores', {})
                    st.markdown(f"- **Overall Score:** `{scores.get('composite', 0):.3f}`")
                    st.markdown(f"- **Doench 2016:** {scores.get('doench_2016', 0):.3f}")
                    st.markdown(f"- **Moreno-Mateos:** {scores.get('moreno_mateos', 0):.3f}")
                    st.markdown(f"- **Xu Score:** {scores.get('xu', 0):.3f}")
                    st.markdown(f"- **GC Score:** {scores.get('gc', 0):.3f}")
                
                # Sequence context
                if grna_data.get('context_sequence'):
                    st.markdown("**Sequence Context (±10bp)**")
                    context = grna_data['context_sequence']
                    # Highlight the gRNA within context
                    st.code(context, language=None)
                
                # Off-target info if available
                if grna_data.get('off_targets'):
                    st.markdown("**Off-Target Analysis**")
                    off_targets = grna_data['off_targets']
                    st.markdown(f"- Off-target sites: {len(off_targets)}")
                
                # Base editing targets if available
                if grna_data.get('base_editing_targets'):
                    st.markdown("**Base Editing Targets**")
                    st.markdown(f"- Edit positions: {grna_data['base_editing_targets']}")
                
                # Copy buttons
                col_copy1, col_copy2, col_copy3 = st.columns(3)
                with col_copy1:
                    st.code(seq, language=None)
                    st.caption("gRNA sequence (copy above)")
                with col_copy2:
                    full_target = f"{seq}{grna_data.get('pam', 'NGG')}"
                    st.code(full_target, language=None)
                    st.caption("Full target + PAM")
                with col_copy3:
                    # Reverse complement for ordering
                    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
                    rev_comp = ''.join(complement.get(b, b) for b in reversed(seq))
                    st.code(rev_comp, language=None)
                    st.caption("Reverse complement")
        
        st.markdown("---")
        
        # Save options
        st.subheader("💾 Save Selected gRNAs")
        col1, col2, col3 = st.columns(3)
        
        with col1:
            gene_folder = st.text_input("Gene folder name", 
                                       value=st.session_state.grna_results.get('metadata', {}).get('gene_symbol', '').lower(),
                                       help="Folder name in data directory",
                                       key="grna_gene_folder")
        
        with col2:
            save_format = st.selectbox("Save format", ["sequence", "fasta", "detailed"], key="grna_save_format")
        
        with col3:
            if st.button("💾 Save Selected gRNAs", key="grna_save_button"):
                save_grnas_to_file(selected_grnas, gene_folder, save_format)


def save_grnas_to_file(grna_sequences, gene_folder, format_type="sequence"):
    """Save selected gRNAs to file."""
    if not gene_folder:
        st.error("Please specify a gene folder name")
        return
    
    # Get data directory - use configured or default to current working directory
    data_dir = st.session_state.get('data_directory', '')
    
    if not data_dir:
        # Use a default output directory in the current working directory
        data_dir = os.path.join(os.getcwd(), "output", "grna_designs")
        st.warning(f"⚠️ No data directory configured. Saving to: `{data_dir}`")
    
    # Create folder structure
    gene_path = os.path.join(data_dir, gene_folder.lower())
    
    try:
        os.makedirs(gene_path, exist_ok=True)
        os.makedirs(os.path.join(gene_path, "input"), exist_ok=True)
    except Exception as e:
        st.error(f"❌ Cannot create directory: {gene_path}")
        st.error(f"Error: {str(e)}")
        st.info("💡 Please set a valid data directory in the sidebar")
        return
    
    # Save gRNA file
    grna_file = os.path.join(gene_path, "grna.txt")
    
    try:
        with open(grna_file, 'w') as f:
            if format_type == "sequence":
                for seq in grna_sequences:
                    f.write(seq + '\n')
            elif format_type == "fasta":
                for i, seq in enumerate(grna_sequences):
                    f.write(f">{gene_folder}_gRNA_{i+1}\n")
                    f.write(seq + '\n')
            elif format_type == "detailed":
                f.write(f"# gRNAs for {gene_folder}\n")
                f.write(f"# Generated by CasPINS - Cas-Primer-Indel Suite\n")
                f.write(f"# Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                for i, seq in enumerate(grna_sequences):
                    f.write(f"# gRNA {i+1}\n")
                    f.write(seq + '\n\n')
        
        st.success(f"✅ Saved {len(grna_sequences)} gRNAs to `{grna_file}`")
        st.info("You can now use these gRNAs in the Primer Design and Indel Analysis tabs")
        
        # Also offer download option
        with open(grna_file, 'r') as f:
            content = f.read()
        st.download_button(
            label="📥 Download gRNA file",
            data=content,
            file_name=f"{gene_folder}_grna.txt",
            mime="text/plain"
        )
        
    except Exception as e:
        st.error(f"Error saving file: {str(e)}")


def primer_design_tab():
    """Primer Design interface - Completely Generic and Independent."""
    st.header("🧪 Primer Design")
    st.markdown("Design PCR primers for any gene from any species - **Completely Independent Tool**")
    
    # Initialize variables
    target_sequence = ""
    grna_sequences = []
    gene_folder = None
    
    # Input method selection
    st.subheader("🎯 Target Selection")
    
    # Auto-import button if data exists in session state from gRNA design
    if st.session_state.get('grna_gene_input') or st.session_state.get('grna_species_expanded'):
        if st.button("🔄 Import settings and gRNAs from gRNA Design", width="stretch", help="Load gene, species, and selected gRNAs from the design tab"):
            if 'grna_gene_input' in st.session_state:
                st.session_state.primer_gene_generic = st.session_state.grna_gene_input
            if 'grna_species_expanded' in st.session_state:
                st.session_state.primer_species_generic = st.session_state.grna_species_expanded
            
            selected_grnas = st.session_state.get('selected_grnas', [])
            if selected_grnas:
                # selected_grnas is a list of plain sequence strings
                grna_strings = [g['sequence'] if isinstance(g, dict) else g for g in selected_grnas]
                st.session_state.primer_grna_input = '\n'.join(grna_strings)
                st.session_state.primer_input_method = "🧬 Paste/Upload Sequence"
            st.rerun()
            
    col1, col2 = st.columns([2, 1])
    
    with col1:
        # Generic gene input (not limited to local files)
        gene_name = st.text_input(
            "Gene Name/Symbol", 
            placeholder="e.g., TP53, GAPDH, BRCA1, Slc6a3, etc.",
            help="Enter any gene name or symbol from any species",
            key="primer_gene_generic"
        )
    
    with col2:
        # Use same comprehensive species list as gRNA design (90+ species)
        primer_species_options = get_species_for_dropdown()
        primer_species_options.append("Other (specify below)")
        
        species = st.selectbox(
            "Species", 
            primer_species_options, 
            index=0,
            help="Select from 90+ supported species (same as gRNA design)",
            key="primer_species_generic"
        )
        
        # Custom species input if "Other" selected
        if species == "Other (specify below)":
            custom_species = st.text_input(
                "Custom Species",
                placeholder="e.g., Macaca fascicularis",
                key="primer_custom_species"
            )
            if custom_species:
                species = custom_species

    # Sequence input methods
    st.subheader("📄 Sequence Input")
    
    # Input method selection
    input_method = st.radio(
        "How would you like to provide the target sequence?",
        ["🧬 Paste/Upload Sequence", "🔍 Search Online Database", "📁 Use Local File"],
        help="Choose your preferred method for sequence input",
        key="primer_input_method"
    )
    
    if input_method == "🧬 Paste/Upload Sequence":
        # Direct sequence input
        col1, col2 = st.columns([2, 1])
        
        with col1:
            target_sequence = st.text_area(
                "Target Sequence (DNA/mRNA)",
                height=200,
                placeholder=">Gene_sequence\nATGGCGACCTTCGACGTGCAGCTGGTGGAGTTC...",
                help="Paste FASTA format or raw sequence. Can be genomic DNA, mRNA, or cDNA",
                key="primer_sequence_input"
            )
        
        with col2:
            # gRNA input - Essential for optimal primer design
            st.markdown("**🧬 gRNA Sequences (Recommended)**")
            grna_input = st.text_area(
                "gRNA Sequences",
                height=100,
                placeholder="ACGTACGTACGTACGTACGT\nGCGTGCGTGCGTGCGTGCGT",
                help="One gRNA per line (20nt each). Essential for positioning primers around the cut site.",
                key="primer_grna_input"
            )
            
            if grna_input:
                grna_sequences = [line.strip().upper() for line in grna_input.split('\n') 
                                if line.strip() and not line.startswith('#')]
                st.success(f"✅ {len(grna_sequences)} gRNA(s) entered")
            else:
                st.warning("⚠️ No gRNA - primers will cover entire sequence")
    
    elif input_method == "🔍 Search Online Database":
        # Online database search - NOW IMPLEMENTED
        st.success("🔗 **Database Integration Active**: Search NCBI, Ensembl for gene sequences")
        
        col1, col2 = st.columns([2, 1])
        
        with col1:
            search_id = st.text_input(
                "Gene ID/Accession", 
                placeholder="e.g., NM_003054.5, ENSG00000107669, DDC",
                help="Enter RefSeq ID (NM_*), Ensembl ID (ENSG*), or gene symbol",
                key="primer_search_id"
            )
        
        with col2:
            search_db = st.selectbox(
                "Database", 
                ["Auto-detect", "NCBI RefSeq", "Ensembl"],
                key="primer_search_db"
            )
        
        # gRNA input - ESSENTIAL for primer design
        st.markdown("#### 🧬 gRNA Sequences (Required for optimal primer design)")
        st.warning("⚠️ **gRNA sequences are essential** for designing primers that flank the cut site correctly.")
        
        grna_input_online = st.text_area(
            "Enter your gRNA sequences",
            height=100,
            placeholder="Enter one gRNA sequence per line (20nt each):\nACGTACGTACGTACGTACGT\nGCGTGCGTGCGTGCGTGCGT",
            help="gRNA sequences determine the cut site location. Primers will be designed to flank this region.",
            key="primer_grna_online"
        )
        
        if grna_input_online:
            grna_sequences = [line.strip().upper() for line in grna_input_online.split('\n') 
                            if line.strip() and not line.startswith('#')]
            st.success(f"✅ {len(grna_sequences)} gRNA sequence(s) entered")
        else:
            st.info("💡 Tip: First design gRNAs in the **gRNA Design** tab, then copy them here.")
        
        # Sequence type selection - IMPORTANT for gRNA compatibility
        st.markdown("#### 📊 Sequence Type")
        seq_type = st.radio(
            "Select sequence type to fetch",
            ["Genomic DNA (for gRNA/CRISPR primers)", "cDNA/mRNA (for expression primers)"],
            index=0,
            help="**Genomic DNA**: Use this when you have gRNA sequences - gRNAs target genomic DNA including introns. "
                 "**cDNA/mRNA**: Use this for primers targeting expressed transcripts only.",
            key="primer_seq_type_online"
        )
        
        use_genomic = "Genomic" in seq_type
        
        if st.button("🔍 Search Database", key="primer_db_search_btn"):
            if search_id:
                with st.spinner(f"Fetching sequence from {search_db}..."):
                    try:
                        import requests
                        
                        fetched_sequence = None
                        sequence_info = {}
                        
                        # Auto-detect database based on ID format
                        if search_db == "Auto-detect":
                            if search_id.startswith("NM_") or search_id.startswith("NR_") or search_id.startswith("NG_"):
                                search_db = "NCBI RefSeq"
                            elif search_id.startswith("ENS"):
                                search_db = "Ensembl"
                            else:
                                # Try gene symbol lookup via Ensembl first
                                search_db = "Ensembl"
                        
                        if search_db == "NCBI RefSeq":
                            # Fetch from NCBI
                            base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
                            efetch_url = f"{base_url}/efetch.fcgi"
                            params = {
                                'db': 'nuccore',
                                'id': search_id,
                                'rettype': 'fasta',
                                'retmode': 'text'
                            }
                            
                            response = requests.get(efetch_url, params=params, timeout=30)
                            
                            if response.status_code == 200 and len(response.text) > 50:
                                lines = response.text.strip().split('\n')
                                header = lines[0] if lines else ''
                                fetched_sequence = ''.join(lines[1:]).upper()
                                sequence_info = {
                                    'source': 'NCBI RefSeq',
                                    'accession': search_id,
                                    'description': header[1:] if header.startswith('>') else header
                                }
                            else:
                                st.error(f"Could not fetch {search_id} from NCBI. Check the accession number.")
                        
                        elif search_db == "Ensembl":
                            # Convert species selection to Ensembl format
                            ensembl_species = convert_species_to_ensembl(species)
                            logger.info(f"Using Ensembl species: {ensembl_species} (from: {species})")
                            st.info(f"🔎 Searching in species: **{ensembl_species}**")
                            
                            url = None
                            
                            # Handle different Ensembl ID types
                            if search_id.upper().startswith("ENSG"):
                                # Gene ID - need to get transcript first
                                st.info(f"🔍 Looking up gene {search_id}...")
                                print(f"[LOG] Fetching gene info for {search_id}")
                                
                                gene_url = f"https://rest.ensembl.org/lookup/id/{search_id}?expand=1"
                                gene_resp = requests.get(gene_url, headers={"Content-Type": "application/json"}, timeout=30)
                                
                                print(f"[LOG] Gene lookup status: {gene_resp.status_code}")
                                
                                if gene_resp.status_code == 200:
                                    gene_data = gene_resp.json()
                                    sequence_info['gene_name'] = gene_data.get('display_name', search_id)
                                    sequence_info['ensembl_id'] = search_id
                                    sequence_info['chromosome'] = gene_data.get('seq_region_name', '')
                                    sequence_info['description'] = gene_data.get('description', '')
                                    
                                    # Check if we should fetch genomic or cDNA
                                    if use_genomic:
                                        # Fetch GENOMIC sequence (includes introns) - needed for gRNA primers
                                        chrom = gene_data.get('seq_region_name', '')
                                        start = gene_data.get('start', 1)
                                        end = gene_data.get('end', start + 1000)
                                        
                                        # Add flanking regions for primer design
                                        flank = 500
                                        seq_start = max(1, start - flank)
                                        seq_end = end + flank
                                        
                                        sequence_info['genomic_start'] = seq_start
                                        sequence_info['genomic_end'] = seq_end
                                        sequence_info['sequence_type'] = 'genomic'
                                        
                                        url = f"https://rest.ensembl.org/sequence/region/{ensembl_species}/{chrom}:{seq_start}..{seq_end}:1"
                                        st.info(f"📋 Fetching genomic sequence: {chrom}:{seq_start}-{seq_end}")
                                        print(f"[LOG] Fetching genomic: {chrom}:{seq_start}..{seq_end}")
                                    else:
                                        # Fetch cDNA sequence (spliced)
                                        transcripts = gene_data.get('Transcript', [])
                                        print(f"[LOG] Found {len(transcripts)} transcripts")
                                        
                                        canonical = next((t for t in transcripts if t.get('is_canonical')), 
                                                        transcripts[0] if transcripts else None)
                                        
                                        if canonical:
                                            transcript_id = canonical.get('id', '')
                                            sequence_info['transcript_id'] = transcript_id
                                            sequence_info['sequence_type'] = 'cdna'
                                            url = f"https://rest.ensembl.org/sequence/id/{transcript_id}?type=cdna"
                                            st.info(f"📋 Found canonical transcript: {transcript_id}")
                                            print(f"[LOG] Using transcript {transcript_id}")
                                        else:
                                            st.error(f"No transcripts found for gene {search_id}")
                                            print(f"[LOG] No transcripts found")
                                else:
                                    error_msg = gene_resp.text[:200] if gene_resp.text else "Unknown error"
                                    st.error(f"Gene {search_id} not found in Ensembl (Status: {gene_resp.status_code})")
                                    print(f"[LOG] Gene lookup failed: {error_msg}")
                            
                            elif search_id.upper().startswith("ENST"):
                                # Transcript ID - direct lookup
                                st.info(f"🔍 Looking up transcript {search_id}...")
                                url = f"https://rest.ensembl.org/sequence/id/{search_id}?type=cdna"
                                sequence_info['transcript_id'] = search_id
                                print(f"[LOG] Direct transcript lookup: {search_id}")
                            
                            else:
                                # Gene symbol - search first
                                st.info(f"🔍 Searching for gene symbol '{search_id}' in {ensembl_species}...")
                                print(f"[LOG] Searching for symbol {search_id} in {ensembl_species}")
                                
                                search_url = f"https://rest.ensembl.org/lookup/symbol/{ensembl_species}/{search_id}"
                                search_resp = requests.get(search_url, headers={"Content-Type": "application/json"}, timeout=30)
                                
                                print(f"[LOG] Symbol search status: {search_resp.status_code}")
                                
                                if search_resp.status_code == 200:
                                    gene_data = search_resp.json()
                                    gene_id = gene_data.get('id', '')
                                    sequence_info['gene_name'] = gene_data.get('display_name', search_id)
                                    sequence_info['ensembl_id'] = gene_id
                                    sequence_info['chromosome'] = gene_data.get('seq_region_name', '')
                                    
                                    st.info(f"📋 Found gene: {gene_id}")
                                    
                                    if use_genomic:
                                        # Fetch GENOMIC sequence for gRNA primers
                                        chrom = gene_data.get('seq_region_name', '')
                                        start = gene_data.get('start', 1)
                                        end = gene_data.get('end', start + 1000)
                                        
                                        flank = 500
                                        seq_start = max(1, start - flank)
                                        seq_end = end + flank
                                        
                                        sequence_info['genomic_start'] = seq_start
                                        sequence_info['genomic_end'] = seq_end
                                        sequence_info['sequence_type'] = 'genomic'
                                        
                                        url = f"https://rest.ensembl.org/sequence/region/{ensembl_species}/{chrom}:{seq_start}..{seq_end}:1"
                                        st.info(f"📋 Fetching genomic sequence: {chrom}:{seq_start}-{seq_end}")
                                        print(f"[LOG] Fetching genomic: {chrom}:{seq_start}..{seq_end}")
                                    else:
                                        # Get transcript info for cDNA
                                        tx_url = f"https://rest.ensembl.org/lookup/id/{gene_id}?expand=1"
                                        tx_resp = requests.get(tx_url, headers={"Content-Type": "application/json"}, timeout=30)
                                        
                                        if tx_resp.status_code == 200:
                                            tx_data = tx_resp.json()
                                            transcripts = tx_data.get('Transcript', [])
                                            canonical = next((t for t in transcripts if t.get('is_canonical')), 
                                                            transcripts[0] if transcripts else None)
                                            if canonical:
                                                transcript_id = canonical.get('id', '')
                                                url = f"https://rest.ensembl.org/sequence/id/{transcript_id}?type=cdna"
                                                sequence_info['transcript_id'] = transcript_id
                                                sequence_info['sequence_type'] = 'cdna'
                                                print(f"[LOG] Using transcript {transcript_id}")
                                            else:
                                                st.error(f"No transcripts found for {search_id}")
                                else:
                                    st.error(f"Gene symbol '{search_id}' not found in Ensembl for {ensembl_species}")
                                    print(f"[LOG] Symbol not found: {search_resp.text[:200]}")
                            
                            # Fetch sequence if we have a URL
                            if url:
                                print(f"[LOG] Fetching sequence from: {url}")
                                seq_resp = requests.get(url, headers={"Content-Type": "text/plain"}, timeout=30)
                                
                                print(f"[LOG] Sequence fetch status: {seq_resp.status_code}")
                                
                                if seq_resp.status_code == 200:
                                    fetched_sequence = seq_resp.text.strip().upper()
                                    sequence_info['source'] = 'Ensembl'
                                    print(f"[LOG] Got sequence: {len(fetched_sequence)} bp")
                                else:
                                    error_detail = seq_resp.text[:200] if seq_resp.text else "Unknown"
                                    st.error(f"Could not fetch sequence (Status: {seq_resp.status_code})")
                                    print(f"[LOG] Sequence fetch failed: {error_detail}")
                        
                        # Store fetched sequence in session state
                        if fetched_sequence and len(fetched_sequence) > 50:
                            st.session_state['fetched_sequence'] = fetched_sequence
                            st.session_state['fetched_sequence_info'] = sequence_info
                            st.success(f"✅ Found sequence: {len(fetched_sequence)} bp from {sequence_info.get('source', 'database')}")
                            
                            # Show sequence info
                            with st.expander("📋 Sequence Information", expanded=True):
                                for key, value in sequence_info.items():
                                    if value and key != 'source':
                                        st.markdown(f"**{key.replace('_', ' ').title()}:** {value}")
                                st.text(f"Sequence (first 200bp): {fetched_sequence[:200]}...")
                        else:
                            st.warning("No valid sequence found. Please try a different ID.")
                            
                    except requests.exceptions.Timeout:
                        st.error("Database request timed out. Please try again.")
                    except Exception as e:
                        st.error(f"Error fetching sequence: {str(e)}")
            else:
                st.warning("Please enter a Gene ID or Accession number")
        
        # Use fetched sequence if available
        if st.session_state.get('fetched_sequence'):
            target_sequence = st.session_state['fetched_sequence']
            st.success(f"✅ Using fetched sequence ({len(target_sequence)} bp)")
            
            # Show sequence preview
            with st.expander("📋 Sequence Preview", expanded=False):
                st.text(f"First 300bp: {target_sequence[:300]}...")
            
            # Remind about gRNA if not entered
            if not grna_sequences:
                st.warning("⚠️ No gRNA sequences entered above. Primers will be designed for the entire sequence.")
    
    else:  # Use Local File
        # Check for local files (backward compatibility)
        if gene_name:
            data_dir = st.session_state.get('data_directory', '')
            gene_folder = os.path.join(data_dir, gene_name.lower()) if data_dir else ""
            
            if os.path.exists(gene_folder):
                # Check for existing files
                grna_file = os.path.join(gene_folder, "grna.txt")
                mrna_file = os.path.join(gene_folder, "mrna.txt")
                
                has_grna = os.path.exists(grna_file)
                has_mrna = os.path.exists(mrna_file)
                
                if has_grna:
                    st.success(f"✓ Found gRNA file for {gene_name}")
                    with open(grna_file, 'r') as f:
                        grna_sequences = [line.strip().upper() for line in f 
                                        if line.strip() and not line.startswith('#')]
                
                if has_mrna:
                    st.success(f"✓ Found mRNA file for {gene_name}")
                    with open(mrna_file, 'r') as f:
                        target_sequence = ''.join(line.strip() for line in f 
                                                if not line.startswith('>')).upper()
            else:
                st.warning(f"No local files found for {gene_name}")
                st.info(f"Set data directory in sidebar, then create folder `{gene_name.lower()}/` with `grna.txt` and `mrna.txt` files")
    
    # Primer design settings
    if target_sequence or (gene_folder and os.path.exists(gene_folder)):
        st.subheader("⚙️ Primer Design Parameters")
        
        # Primer set selection
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("#### 🧬 PCR I - Genomic DNA Primers")
            st.info("Large flanking regions (~500bp) for genomic DNA amplification")
            pcr1_enabled = st.checkbox("Design PCR I primers", value=True, key="primer_pcr1_generic")
            
        with col2:
            st.markdown("#### 🔬 PCR II - Sequencing Primers")
            st.info("Optimized flanking (~200bp) for Sanger sequencing")
            pcr2_enabled = st.checkbox("Design PCR II primers", value=True, key="primer_pcr2_generic")
        
        # Advanced primer parameters
        with st.expander("🔧 Advanced Primer Parameters"):
            col1, col2, col3 = st.columns(3)
            
            with col1:
                st.markdown("**Temperature (°C)**")
                min_tm = st.number_input("Min Tm", min_value=45.0, max_value=75.0, value=57.0, step=0.5, key="primer_min_tm_generic")
                max_tm = st.number_input("Max Tm", min_value=45.0, max_value=75.0, value=63.0, step=0.5, key="primer_max_tm_generic")
                opt_tm = st.number_input("Optimal Tm", min_value=45.0, max_value=75.0, value=60.0, step=0.5, key="primer_opt_tm_generic")
            
            with col2:
                st.markdown("**GC Content (%)**")
                min_gc = st.number_input("Min GC%", min_value=20, max_value=80, value=40, key="primer_min_gc_generic")
                max_gc = st.number_input("Max GC%", min_value=20, max_value=80, value=60, key="primer_max_gc_generic")
                opt_gc = st.number_input("Optimal GC%", min_value=20, max_value=80, value=50, key="primer_opt_gc_generic")
            
            with col3:
                st.markdown("**Primer Length**")
                primer_len_range = st.slider("Length range", 18, 30, (20, 25), key="primer_len_generic")
                max_poly_x = st.number_input("Max poly-X", min_value=3, max_value=8, value=5, key="primer_max_poly_generic")
                
        # Additional options
        with st.expander("🎯 Target Region Options"):
            if grna_sequences:
                st.success(f"Found {len(grna_sequences)} gRNA sequence(s)")
                for i, grna in enumerate(grna_sequences):
                    st.text(f"gRNA {i+1}: {grna}")
            else:
                st.info("No gRNA sequences provided - will design primers for entire sequence")
            
            # Manual target region specification
            manual_target = st.checkbox("Specify custom target region", key="primer_manual_target")
            if manual_target:
                col1, col2 = st.columns(2)
                with col1:
                    target_start = st.number_input("Target start position", min_value=1, value=1, key="primer_target_start")
                with col2:
                    target_end = st.number_input("Target end position", min_value=1, value=100, key="primer_target_end")
        
        # Design button
        if st.button("🚀 Design Primers", type="primary", width="stretch", key="primer_design_generic"):
            if not target_sequence:
                st.error("Please provide a target sequence")
                return
            
            design_primers_generic(
                gene_name, species, target_sequence, grna_sequences,
                pcr1_enabled, pcr2_enabled, 
                min_tm, max_tm, opt_tm, min_gc, max_gc, opt_gc,
                primer_len_range, max_poly_x
            )
    
    elif gene_name:
        st.info("Please provide a target sequence to design primers")
        st.markdown("""
        **Options:**
        1. **Paste sequence** directly in the text area above
        2. **Upload a FASTA file** with your sequence
        3. **Use local files** if you have them in your configured data directory
        """)


def design_primers_generic(gene_name, species, target_sequence, grna_sequences, 
                          pcr1_enabled, pcr2_enabled, min_tm, max_tm, opt_tm, 
                          min_gc, max_gc, opt_gc, primer_len_range, max_poly_x):
    """Run generic primer design for any gene/species using the enhanced system."""
    
    with st.spinner(f"Designing primers for {gene_name} ({species})..."):
        try:
            # Clean up sequence
            if target_sequence.startswith('>'):
                # Remove FASTA header
                lines = target_sequence.split('\n')
                target_sequence = ''.join(line.strip() for line in lines[1:])
            
            target_sequence = target_sequence.upper().replace(' ', '').replace('\n', '')
            
            # Validate sequence
            valid_bases = set('ATCGRYSWKMBDHVN')
            if not all(base in valid_bases for base in target_sequence):
                st.error("Invalid characters in sequence. Please use standard DNA bases.")
                return
            
            if len(target_sequence) < 100:
                st.error("Sequence too short. Please provide at least 100 bases.")
                return
            
            st.markdown(f"### 🧬 Enhanced Primer Design Results for {gene_name}")
            st.markdown(f"**Species**: {species}")
            st.markdown(f"**Target sequence length**: {len(target_sequence)} bp")
            
            # Find gRNA positions if provided
            grna_positions = []
            if grna_sequences:
                st.markdown(f"**Found {len(grna_sequences)} gRNA sequence(s)**")
                
                for i, grna in enumerate(grna_sequences):
                    pos, strand, cut_site = find_grna_in_sequence(target_sequence, grna)
                    if cut_site is not None:
                        grna_positions.append((pos, strand, cut_site))
                        st.text(f"gRNA {i+1}: {grna}")
                        st.text(f"  Position: {pos}, Strand: {strand}, Cut site: {cut_site}")
                    else:
                        st.warning(f"gRNA {i+1}: {grna} - NOT FOUND in target sequence!")
            
            # Convert species name to simple format
            species_simple = species.split('(')[0].strip().lower()
            if 'human' in species_simple or 'homo' in species_simple:
                species_simple = 'human'
            elif 'mouse' in species_simple or 'mus' in species_simple:
                species_simple = 'mouse'
            elif 'rat' in species_simple or 'rattus' in species_simple:
                species_simple = 'rat'
            else:
                species_simple = 'human'  # Default fallback
            
            # Use enhanced primer design system
            if grna_positions:
                # Generate enhanced primer report
                primer_report = generate_enhanced_primer_report(
                    grna_positions, target_sequence, gene_name, species_simple
                )
                
                # Display the full report
                st.markdown("### 📋 Comprehensive Primer Design Report")
                st.text(primer_report)
                
                # Also create a more interactive display
                display_enhanced_primer_results(grna_positions, target_sequence, gene_name, species_simple)
                
                # Download options
                st.markdown("### 📥 Download Options")
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    st.download_button(
                        label="📄 Download Full Report",
                        data=primer_report,
                        file_name=f"primer_report_{gene_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt",
                        mime="text/plain",
                        key="download_full_report"
                    )
                
                with col2:
                    # Generate ready-to-order format
                    order_format = generate_order_format(primer_report, gene_name)
                    st.download_button(
                        label="🛒 Ready-to-Order Format",
                        data=order_format,
                        file_name=f"primers_to_order_{gene_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt",
                        mime="text/plain",
                        key="download_order_format"
                    )
                
                with col3:
                    # Generate JSON format
                    json_data = generate_json_format(grna_positions, target_sequence, gene_name, species_simple)
                    st.download_button(
                        label="📊 JSON Data",
                        data=json.dumps(json_data, indent=2),
                        file_name=f"primer_data_{gene_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                        mime="application/json",
                        key="download_json_data"
                    )
            else:
                st.warning("No gRNA sequences provided or found. Using fallback primer design.")
                # Fallback to basic primer design for middle region
                seq_len = len(target_sequence)
                middle_pos = seq_len // 2
                fake_grna_positions = [(middle_pos - 10, 1, middle_pos)]
                
                primer_report = generate_enhanced_primer_report(
                    fake_grna_positions, target_sequence, gene_name, species_simple
                )
                
                st.text(primer_report)
                
                st.download_button(
                    label="📄 Download Report",
                    data=primer_report,
                    file_name=f"primer_report_{gene_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt",
                    mime="text/plain",
                    key="download_fallback_report"
                )
            
        except Exception as e:
            st.error(f"Error during primer design: {str(e)}")
            st.error("Please check your sequence and try again.")
            # Show detailed error for debugging
            with st.expander("Debug Information"):
                st.code(str(e))


def display_enhanced_primer_results(grna_positions, target_sequence, gene_name, species):
    """Display enhanced primer results in an interactive format."""
    
    st.markdown("### 🎯 Interactive Primer Results")
    
    # Initialize the enhanced designer
    designer = PrimerDesigner(gene_name, species)
    
    # Define target regions
    cut_sites = [pos[2] for pos in grna_positions]
    min_cut = min(cut_sites)
    max_cut = max(cut_sites)
    target_regions = [(min_cut - 50, max_cut + 50)]
    
    # Design comprehensive primers
    results = designer.design_comprehensive_primers(target_sequence, target_regions, 
                                                   {'positions': grna_positions})
    
    # Display PCR I primers
    if results.get('pcr1_primers'):
        st.markdown("#### 🧬 PCR I - Genomic DNA Primers")
        pcr1_data = []
        for primer in results['pcr1_primers'][:3]:  # Top 3
            pcr1_data.append({
                'Set': primer['rank'],
                'Forward': primer['forward_seq'],
                'Reverse': primer['reverse_seq'],
                'F_Tm': f"{primer['forward_tm']}°C",
                'R_Tm': f"{primer['reverse_tm']}°C",
                'Product': f"{primer['product_size']} bp",
                'Quality': f"{primer.get('penalty', 'N/A')}"
            })
        
        pcr1_df = pd.DataFrame(pcr1_data)
        st.dataframe(pcr1_df, width="stretch", hide_index=True)
    
    # Display PCR II primers
    if results.get('pcr2_primers'):
        st.markdown("#### 🔬 PCR II - Sequencing Primers")
        pcr2_data = []
        for primer in results['pcr2_primers'][:3]:  # Top 3
            pcr2_data.append({
                'Set': primer['rank'],
                'Forward': primer['forward_seq'],
                'Reverse': primer['reverse_seq'],
                'F_Tm': f"{primer['forward_tm']}°C",
                'R_Tm': f"{primer['reverse_tm']}°C",
                'Product': f"{primer['product_size']} bp",
                'Quality': f"{primer.get('penalty', 'N/A')}"
            })
        
        pcr2_df = pd.DataFrame(pcr2_data)
        st.dataframe(pcr2_df, width="stretch", hide_index=True)
    
    # Display quality metrics
    if results.get('quality_metrics'):
        st.markdown("#### 📊 Quality Metrics")
        metrics = results['quality_metrics']
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Total Primers", metrics.get('total_primers', 0))
        with col2:
            st.metric("Avg Quality Score", f"{metrics.get('avg_penalty', 0):.2f}")
        with col3:
            tm_range = metrics.get('tm_range', {})
            st.metric("Tm Range", f"{tm_range.get('min', 0):.1f}-{tm_range.get('max', 0):.1f}°C")
        
        # Show potential issues
        if metrics.get('potential_issues'):
            st.warning("⚠️ Potential Issues Detected:")
            for issue in metrics['potential_issues']:
                st.text(f"• {issue}")
    
    # Display recommendations
    if results.get('recommendations'):
        st.markdown("#### 💡 Recommendations")
        for rec in results['recommendations']:
            if rec.startswith('✓'):
                st.success(rec)
            elif rec.startswith('⚠'):
                st.warning(rec)
            else:
                st.info(rec)


def generate_order_format(primer_report, gene_name):
    """Generate ready-to-order primer format from the report."""
    import re
    
    order_lines = []
    order_lines.append(f"# Primers for {gene_name} - Ready to Order Format")
    order_lines.append(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    order_lines.append(f"# Copy and paste these sequences into your primer ordering form\n")
    
    # Extract PCR I primers
    pcr1_matches = re.findall(r'PCR I - Set (\d+).*?Forward primer: 5\'-([ATCG]+)-3\'.*?Reverse primer: 5\'-([ATCG]+)-3\'', 
                             primer_report, re.DOTALL)
    
    if pcr1_matches:
        order_lines.append("# PCR I Primers (Genomic DNA Amplification)")
        for i, (set_num, fwd, rev) in enumerate(pcr1_matches[:2]):  # Top 2 sets
            order_lines.append(f"{gene_name}_PCR1_F{set_num}\t{fwd}")
            order_lines.append(f"{gene_name}_PCR1_R{set_num}\t{rev}")
        order_lines.append("")
    
    # Extract PCR II primers
    pcr2_matches = re.findall(r'PCR II - Set (\d+).*?Forward primer: 5\'-([ATCG]+)-3\'.*?Reverse primer: 5\'-([ATCG]+)-3\'', 
                             primer_report, re.DOTALL)
    
    if pcr2_matches:
        order_lines.append("# PCR II Primers (Sequencing)")
        for i, (set_num, fwd, rev) in enumerate(pcr2_matches[:3]):  # Top 3 sets
            order_lines.append(f"{gene_name}_PCR2_F{set_num}\t{fwd}")
            order_lines.append(f"{gene_name}_PCR2_R{set_num}\t{rev}")
        order_lines.append("")
    
    # Add ordering notes
    order_lines.extend([
        "\n# ORDERING NOTES:",
        "# - Standard desalting is sufficient for most applications",
        "# - Consider ordering at 25nmol scale",
        "# - For sequencing primers (PCR II), consider PAGE purification",
        "# - Store primers at -20°C upon arrival"
    ])
    
    return "\n".join(order_lines)


def generate_json_format(grna_positions, target_sequence, gene_name, species):
    """Generate JSON format data for machine processing."""
    
    # Initialize the enhanced designer
    designer = PrimerDesigner(gene_name, species)
    
    # Define target regions
    cut_sites = [pos[2] for pos in grna_positions]
    min_cut = min(cut_sites)
    max_cut = max(cut_sites)
    target_regions = [(min_cut - 50, max_cut + 50)]
    
    # Design comprehensive primers
    results = designer.design_comprehensive_primers(target_sequence, target_regions, 
                                                   {'positions': grna_positions})
    
    # Create JSON structure
    json_data = {
        'gene_name': gene_name,
        'species': species,
        'timestamp': datetime.now().isoformat(),
        'target_sequence_length': len(target_sequence),
        'grna_positions': [
            {
                'position': pos[0],
                'strand': pos[1],
                'cut_site': pos[2]
            }
            for pos in grna_positions
        ],
        'pcr1_primers': results.get('pcr1_primers', []),
        'pcr2_primers': results.get('pcr2_primers', []),
        'quality_metrics': results.get('quality_metrics', {}),
        'recommendations': results.get('recommendations', []),
        'database_info': results.get('database_info', {})
    }
    
    return json_data


def design_primers_with_primer3_generic(target_sequence, target_regions, gene_name,
                                       pcr1_enabled, pcr2_enabled, min_tm, max_tm, opt_tm,
                                       min_gc, max_gc, opt_gc, primer_len_range, max_poly_x):
    """Legacy function - now redirects to enhanced system."""
    
    # This function is kept for backward compatibility
    # but now uses the enhanced primer design system
    
    # Create fake gRNA positions for the target regions
    fake_grna_positions = []
    for start, end in target_regions:
        middle = (start + end) // 2
        fake_grna_positions.append((middle - 10, 1, middle))
    
    # Use enhanced system
    species = 'human'  # Default
    primer_report = generate_enhanced_primer_report(
        fake_grna_positions, target_sequence, gene_name, species
    )
    
    # Parse the report to extract primer data (simplified)
    primer_results = {'pcr1': [], 'pcr2': []}
    
    # This is a simplified parser - the enhanced system provides much more detail
    # but we maintain compatibility with the old interface
    
    return primer_results


def display_primer_results(primer_results, gene_name, species):
    """Display primer design results in a user-friendly format."""
    
    st.markdown("### 🎯 Primer Design Results")
    
    # PCR I Results
    if primer_results['pcr1']:
        st.markdown("#### 🧬 PCR I - Genomic DNA Primers")
        st.info("Designed for genomic DNA amplification (larger products, spans introns)")
        
        for i, result in enumerate(primer_results['pcr1']):
            if result.get('PRIMER_PAIR_NUM_RETURNED', 0) > 0:
                with st.expander(f"PCR I Region {i+1} - {result['PRIMER_PAIR_NUM_RETURNED']} primer pairs"):
                    for j in range(result['PRIMER_PAIR_NUM_RETURNED']):
                        left_seq = result[f'PRIMER_LEFT_{j}_SEQUENCE']
                        right_seq = result[f'PRIMER_RIGHT_{j}_SEQUENCE']
                        left_tm = result[f'PRIMER_LEFT_{j}_TM']
                        right_tm = result[f'PRIMER_RIGHT_{j}_TM']
                        product_size = result[f'PRIMER_PAIR_{j}_PRODUCT_SIZE']
                        
                        st.markdown(f"**Pair {j+1}:**")
                        st.code(f"Forward:  5'-{left_seq}-3'  (Tm: {left_tm:.1f}°C)")
                        st.code(f"Reverse:  5'-{right_seq}-3'  (Tm: {right_tm:.1f}°C)")
                        st.text(f"Product size: {product_size} bp")
                        st.markdown("---")
            else:
                st.warning(f"No primers found for PCR I Region {i+1}")
    
    # PCR II Results
    if primer_results['pcr2']:
        st.markdown("#### 🔬 PCR II - Sequencing Primers")
        st.info("Optimized for Sanger sequencing (smaller products, higher specificity)")
        
        for i, result in enumerate(primer_results['pcr2']):
            if result.get('PRIMER_PAIR_NUM_RETURNED', 0) > 0:
                with st.expander(f"PCR II Region {i+1} - {result['PRIMER_PAIR_NUM_RETURNED']} primer pairs"):
                    for j in range(result['PRIMER_PAIR_NUM_RETURNED']):
                        left_seq = result[f'PRIMER_LEFT_{j}_SEQUENCE']
                        right_seq = result[f'PRIMER_RIGHT_{j}_SEQUENCE']
                        left_tm = result[f'PRIMER_LEFT_{j}_TM']
                        right_tm = result[f'PRIMER_RIGHT_{j}_TM']
                        product_size = result[f'PRIMER_PAIR_{j}_PRODUCT_SIZE']
                        
                        st.markdown(f"**Pair {j+1}:**")
                        st.code(f"Forward:  5'-{left_seq}-3'  (Tm: {left_tm:.1f}°C)")
                        st.code(f"Reverse:  5'-{right_seq}-3'  (Tm: {right_tm:.1f}°C)")
                        st.text(f"Product size: {product_size} bp")
                        st.markdown("---")
            else:
                st.warning(f"No primers found for PCR II Region {i+1}")
    
    # Generate downloadable results
    generate_downloadable_primer_results(primer_results, gene_name, species)


def generate_downloadable_primer_results(primer_results, gene_name, species):
    """Generate downloadable primer results file."""
    
    from datetime import datetime
    
    # Generate text report
    report_lines = []
    report_lines.append(f"Primer Design Results for {gene_name} ({species})")
    report_lines.append("=" * 80)
    report_lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append(f"Tool: CasPINS - Cas-Primer-Indel Suite - Generic Primer Designer")
    report_lines.append("")
    
    # PCR I Results
    if primer_results['pcr1']:
        report_lines.append("PCR I - GENOMIC DNA PRIMERS")
        report_lines.append("-" * 40)
        report_lines.append("Purpose: Amplification from genomic DNA (includes introns)")
        report_lines.append("Expected product size: 800-2000 bp")
        report_lines.append("")
        
        for i, result in enumerate(primer_results['pcr1']):
            if result.get('PRIMER_PAIR_NUM_RETURNED', 0) > 0:
                report_lines.append(f"Region {i+1}:")
                for j in range(result['PRIMER_PAIR_NUM_RETURNED']):
                    left_seq = result[f'PRIMER_LEFT_{j}_SEQUENCE']
                    right_seq = result[f'PRIMER_RIGHT_{j}_SEQUENCE']
                    left_tm = result[f'PRIMER_LEFT_{j}_TM']
                    right_tm = result[f'PRIMER_RIGHT_{j}_TM']
                    product_size = result[f'PRIMER_PAIR_{j}_PRODUCT_SIZE']
                    
                    report_lines.append(f"  Pair {j+1}:")
                    report_lines.append(f"    Forward:  5'-{left_seq}-3'  (Tm: {left_tm:.1f}°C)")
                    report_lines.append(f"    Reverse:  5'-{right_seq}-3'  (Tm: {right_tm:.1f}°C)")
                    report_lines.append(f"    Product:  {product_size} bp")
                    report_lines.append("")
    
    # PCR II Results
    if primer_results['pcr2']:
        report_lines.append("PCR II - SEQUENCING PRIMERS")
        report_lines.append("-" * 40)
        report_lines.append("Purpose: Optimized for Sanger sequencing")
        report_lines.append("Expected product size: 300-600 bp")
        report_lines.append("")
        
        for i, result in enumerate(primer_results['pcr2']):
            if result.get('PRIMER_PAIR_NUM_RETURNED', 0) > 0:
                report_lines.append(f"Region {i+1}:")
                for j in range(result['PRIMER_PAIR_NUM_RETURNED']):
                    left_seq = result[f'PRIMER_LEFT_{j}_SEQUENCE']
                    right_seq = result[f'PRIMER_RIGHT_{j}_SEQUENCE']
                    left_tm = result[f'PRIMER_LEFT_{j}_TM']
                    right_tm = result[f'PRIMER_RIGHT_{j}_TM']
                    product_size = result[f'PRIMER_PAIR_{j}_PRODUCT_SIZE']
                    
                    report_lines.append(f"  Pair {j+1}:")
                    report_lines.append(f"    Forward:  5'-{left_seq}-3'  (Tm: {left_tm:.1f}°C)")
                    report_lines.append(f"    Reverse:  5'-{right_seq}-3'  (Tm: {right_tm:.1f}°C)")
                    report_lines.append(f"    Product:  {product_size} bp")
                    report_lines.append("")
    
    report_lines.append("NOTES:")
    report_lines.append("- PCR I primers are designed for genomic DNA and may span introns")
    report_lines.append("- PCR II primers are optimized for cDNA/mRNA sequencing")
    report_lines.append("- Validate primers with your specific PCR conditions")
    report_lines.append("- Consider primer specificity for your target organism")
    
    report_text = "\n".join(report_lines)
    
    # Download button
    st.download_button(
        label="📥 Download Primer Results",
        data=report_text,
        file_name=f"primers_{gene_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt",
        mime="text/plain",
        key="download_primer_results"
    )
    
    st.success("✅ Primer design complete! Use the download button above to save results.")
    st.info("💡 **Next Steps**: Test primers with your specific PCR conditions and validate specificity.")


def indel_analysis_tab():
    """Indel Analysis interface."""
    st.header("📊 Indel Analysis")
    st.markdown("Analyze CRISPR editing outcomes from Sanger sequencing data")
    
    # Analysis mode
    analysis_mode = st.radio("Analysis Mode",
                            ["Single Sample", "Batch Analysis"],
                            horizontal=True, key="indel_analysis_mode")
    
    if analysis_mode == "Single Sample":
        single_sample_analysis()
    else:
        batch_analysis()


def single_sample_analysis():
    """Single sample indel analysis."""
    st.subheader("Single Sample Analysis")
    
    col1, col2 = st.columns(2)
    
    with col1:
        # Gene name — free text, no folder required
        gene_name = st.text_input(
            "Gene Name",
            placeholder="e.g., DDC, TP53, SLC18A2",
            help="Enter gene name. If a matching folder exists in your data directory, gRNA and mRNA will be loaded automatically.",
            key="indel_gene_input"
        )
        
        # File upload
        st.markdown("#### Upload AB1 Files")
        control_file = st.file_uploader("Control AB1 file (unedited)", type=['ab1'], key="indel_control_file")
        edited_file = st.file_uploader("Edited AB1 file", type=['ab1'], key="indel_edited_file")
        
        # gRNA input — auto-filled from folder if available, manual entry otherwise
        st.markdown("#### gRNA Sequence (optional)")
        
        # Try to pre-fill from gene folder
        data_dir = st.session_state.get('data_directory', '')
        prefill_grna = ""
        grna_source = "manual"
        if gene_name and data_dir:
            grna_file = os.path.join(data_dir, gene_name.lower(), "grna.txt")
            if os.path.exists(grna_file):
                with open(grna_file, 'r') as f:
                    loaded = [ln.strip().upper() for ln in f if ln.strip()]
                if loaded:
                    prefill_grna = '\n'.join(loaded)
                    grna_source = "folder"
        
        grna_input = st.text_area(
            "gRNA sequence(s)",
            value=prefill_grna,
            placeholder="ACGTACGTACGTACGTACGT\n(one per line, 20nt each)",
            height=80,
            help="Used to locate the expected cut site in the trace. If left blank, CasPINS will analyze the full trace without cut-site anchoring.",
            key="indel_grna_manual"
        )
        if grna_source == "folder" and prefill_grna:
            st.caption(f"✅ Auto-loaded from `{gene_name.lower()}/grna.txt`")
        elif grna_input:
            st.caption("✏️ Using manually entered gRNA(s)")
        else:
            st.caption("ℹ️ No gRNA provided — analysis will run without cut-site anchoring")
    
    with col2:
        # Show what we found in the gene folder (informational only)
        if gene_name and data_dir:
            gene_folder = os.path.join(data_dir, gene_name.lower())
            if os.path.exists(gene_folder):
                mrna_path = os.path.join(gene_folder, "mrna.txt")
                mrna_exists = os.path.exists(mrna_path)
                st.info(
                    f"📁 Gene folder found: `{gene_name.lower()}/`\n\n"
                    f"• gRNA file: {'✅' if prefill_grna else '❌'}\n"
                    f"• mRNA file: {'✅' if mrna_exists else '❌ (optional)'}"
                )
            else:
                st.info("ℹ️ No gene folder found — analysis will use uploaded files and entered gRNA only.")
        
        st.markdown("#### Analysis Settings")
        r_squared_correction = st.checkbox(
            "Conservative mode (R² correction)",
            value=False,
            help="Multiply editing estimate by R² goodness-of-fit, producing conservative estimates matching TIDE/ICE behavior on noisy data. Default OFF — CasPINS' combined-channel NNLS extracts more signal from noisy traces by design; low R² is reported as a confidence warning rather than suppressing the estimate."
        )
    
    # Analysis button
    if st.button("🔬 Analyze Indels", type="primary", width="stretch", key="indel_analyze_button"):
        if not control_file or not edited_file:
            st.error("Please upload both control and edited AB1 files.")
            return
        if not gene_name:
            gene_name = "unknown"
        
        grna_sequences = [ln.strip().upper() for ln in grna_input.split('\n') if ln.strip()] if grna_input else []
        run_single_sample_analysis(
            control_file, edited_file, gene_name,
            grna_sequences=grna_sequences,
            r_squared_correction=r_squared_correction
        )


def run_single_sample_analysis(control_file, edited_file, gene_name, grna_sequences=None, r_squared_correction=True):
    """Run indel analysis on a single sample with TIDE-style visualization."""
    with st.spinner("Analyzing indels..."):
        try:
            # Save uploaded files temporarily
            with tempfile.NamedTemporaryFile(delete=False, suffix='.ab1') as tmp_control:
                tmp_control.write(control_file.read())
                control_path = tmp_control.name
            
            with tempfile.NamedTemporaryFile(delete=False, suffix='.ab1') as tmp_edited:
                tmp_edited.write(edited_file.read())
                edited_path = tmp_edited.name
            
            # Try to load gRNA from folder if not supplied from UI
            if grna_sequences is None:
                grna_sequences = []
            data_dir = st.session_state.get('data_directory', '')
            gene_folder = os.path.join(data_dir, gene_name.lower()) if data_dir else ""
            
            if not grna_sequences and gene_folder and os.path.exists(gene_folder):
                grna_path = os.path.join(gene_folder, "grna.txt")
                if os.path.exists(grna_path):
                    with open(grna_path, 'r') as f:
                        grna_sequences = [ln.strip().upper() for ln in f if ln.strip()]
            
            # Try to load mRNA from folder (optional — used only for cut-site location)
            mrna_seq = None
            if gene_folder and os.path.exists(gene_folder):
                mrna_path = os.path.join(gene_folder, "mrna.txt")
                if os.path.exists(mrna_path):
                    with open(mrna_path, 'r') as f:
                        mrna_seq = ''.join(ln.strip() for ln in f).upper()
            
            # Parse AB1 files
            control_seq, control_traces = parse_ab1(control_path)
            edited_seq, edited_traces = parse_ab1(edited_path)
            
            # Find cut sites (requires both gRNA sequences and mRNA reference)
            expected_cut_site = None
            if grna_sequences and mrna_seq:
                cut_sites = []
                for grna in grna_sequences:
                    pos, strand, cut_site = find_grna_in_sequence(mrna_seq, grna)
                    if cut_site is not None:
                        cut_sites.append(cut_site)
                expected_cut_site = cut_sites[0] if cut_sites else None
            elif grna_sequences and not mrna_seq:
                st.info("ℹ️ No mRNA reference file found — cut site will not be anchored. Analysis runs on full trace.")
            
            # Perform indel analysis
            try:
                efficiency_data = decompose_traces_indel_analysis(
                    control_traces, edited_traces, expected_cut_site, r_squared_correction=r_squared_correction
                )
                efficiency_data['method'] = 'Trace Decomposition'
            except:
                efficiency_data = calculate_editing_efficiency_fallback(
                    control_seq, control_traces, edited_seq, edited_traces, expected_cut_site
                )
                efficiency_data['method'] = 'Fallback'
            
            # Display results - TIDE STYLE IN GUI using native Streamlit components
            st.success("✅ Analysis Complete!")
            
            # Show R² quality warning when model fit is poor — inform without suppressing
            r_sq = efficiency_data.get('quality_score', 0) / 100.0
            if r_sq < 0.3:
                st.warning(
                    f"⚠️ Low model fit (R² = {r_sq:.3f}). This may reflect inherently noisy Sanger traces "
                    f"(e.g., multiple large indels, poor base-calling). CasPINS reports the full signal extracted "
                    f"from all trace channels. Enable 'Conservative mode' in Analysis Settings for TIDE/ICE-comparable "
                    f"conservative estimates, or consider NGS-based quantification (e.g., CRISPResso2) for samples "
                    f"with R² < 0.3."
                )
            elif r_sq < 0.5:
                st.info(f"ℹ️ Moderate model fit (R² = {r_sq:.3f}). Results are usable but treat efficiency estimate with moderate confidence.")
            
            # Display results using native Streamlit components
            display_results_streamlit(st, efficiency_data, f"{gene_name.upper()} - Single Sample")
            
            # Display indel spectrum using native Streamlit chart
            display_indel_spectrum_streamlit(st, efficiency_data)
            
            # Show confidence level
            confidence = efficiency_data.get('confidence', 'UNKNOWN')
            if 'Fallback' in efficiency_data.get('method', ''):
                st.error(f"🔴 **Confidence: {confidence}** - Fallback method used (signal quality issues)")
            elif confidence == 'HIGH':
                st.success(f"🟢 **Confidence: {confidence}** - Reliable NNLS trace decomposition")
            elif confidence == 'MEDIUM':
                st.warning(f"🟡 **Confidence: {confidence}** - Moderate quality trace decomposition")
            else:
                st.info(f"🔵 **Confidence: {confidence}** - Check trace quality")
            
            # Detailed metrics
            st.markdown("### 📋 Detailed Metrics")
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                st.metric("Editing Efficiency", 
                         f"{efficiency_data.get('editing_efficiency', 0):.1f}%")
            with col2:
                st.metric("Wild-type", 
                         f"{efficiency_data.get('wt_fraction', 0):.1f}%")
            with col3:
                dominant = efficiency_data.get('dominant_indel_size', 'N/A')
                st.metric("Dominant Indel", 
                         f"{dominant} bp" if dominant not in ['N/A', 'Unknown'] else 'N/A')
            with col4:
                st.metric("R² Quality", 
                         f"{efficiency_data.get('quality_score', 0):.1f}%")
            
            # Generate comprehensive analysis plot
            st.markdown("### 🖼️ Indel Analysis Visualization")
            
            # Determine a safe output directory regardless of whether gene_folder exists
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            if gene_folder and os.path.exists(gene_folder):
                output_dir = ensure_output_dir(gene_folder)
            else:
                # Fall back to a session-scoped temp output folder
                fallback_out = os.path.join(
                    st.session_state.get('data_directory', tempfile.gettempdir()),
                    "caspins_output"
                )
                os.makedirs(fallback_out, exist_ok=True)
                output_dir = fallback_out
            
            # Create TIDE-style figure
            try:
                tide_fig = create_tide_style_figure(
                    control_seq, edited_seq,
                    control_traces, edited_traces,
                    expected_cut_site or 200,
                    efficiency_data,
                    grna_sequences[0] if grna_sequences else "N/A",
                    sample_name=f"{gene_name.upper()} - Single Analysis"
                )
                
                # Save and display
                tide_plot_path = os.path.join(output_dir, f"tide_analysis_{gene_name}_{timestamp}.png")
                tide_fig.savefig(tide_plot_path, dpi=150, bbox_inches='tight', facecolor='white')
                
                import matplotlib.pyplot as plt
                plt.close(tide_fig)
                
                st.image(tide_plot_path, caption="Comprehensive Indel Analysis (Sequence Alignment + Chromatogram + Indel Spectrum)")
                
                # Download button for TIDE plot
                with open(tide_plot_path, 'rb') as f:
                    st.download_button(
                        label="📥 Download Analysis Plot",
                        data=f.read(),
                        file_name=f"tide_analysis_{gene_name}_{timestamp}.png",
                        mime="image/png",
                        key="download_tide_plot"
                    )
            except Exception as e:
                st.warning(f"Could not generate TIDE-style plot: {e}")
            
            # Also generate standard plot (optional secondary visualization)
            try:
                plot_path = plot_indel_analysis_multi(
                    control_path, edited_path, output_dir, gene_name,
                    "single_analysis", efficiency_data, grna_sequences, timestamp
                )
                if plot_path and os.path.exists(plot_path):
                    with st.expander("View Standard Analysis Plot"):
                        st.image(plot_path, caption="Standard Indel Analysis Results")
                        with open(plot_path, 'rb') as f:
                            st.download_button(
                                label="📥 Download Standard Plot",
                                data=f.read(),
                                file_name=f"indel_analysis_{gene_name}_{timestamp}.png",
                                mime="image/png",
                                key="download_standard_plot"
                            )
            except Exception:
                pass  # Secondary plot is optional; main TIDE-style plot already shown
            
            # Cleanup temp files
            os.unlink(control_path)
            os.unlink(edited_path)
            
        except Exception as e:
            st.error(f"Error during analysis: {str(e)}")
            import traceback
            with st.expander("Error Details"):
                st.code(traceback.format_exc())


def batch_analysis():
    """Batch analysis interface."""
    st.subheader("Batch Analysis")
    st.info("Process multiple edited samples against a single control")
    
    # Gene selection for batch - use configured data directory
    gene_folders_with_input = []
    data_dir = st.session_state.get('data_directory', '')
    
    if not data_dir or not os.path.isdir(data_dir):
        st.warning("⚠️ Please set a valid data directory in the sidebar settings")
        st.info("💡 The data directory should contain gene folders with input/ subdirectories")
        return
    
    for gene_dir in os.listdir(data_dir):
        gene_path = os.path.join(data_dir, gene_dir)
        input_path = os.path.join(gene_path, "input")
        if os.path.isdir(gene_path) and os.path.exists(input_path):
            # Check for control and edited files
            control_exists = os.path.exists(os.path.join(input_path, "control.ab1"))
            edited_files = len([f for f in os.listdir(input_path) 
                              if f.startswith("edited") and f.endswith(".ab1")])
            grna_exists = os.path.exists(os.path.join(gene_path, "grna.txt"))
            mrna_exists = os.path.exists(os.path.join(gene_path, "mrna.txt"))
            
            if control_exists and edited_files > 0:
                gene_folders_with_input.append({
                        'gene': gene_dir.upper(),
                        'edited_count': edited_files,
                        'ready': grna_exists and mrna_exists,
                        'missing': [] if (grna_exists and mrna_exists) else 
                                  (['grna.txt'] if not grna_exists else []) + 
                                  (['mrna.txt'] if not mrna_exists else [])
                    })
    
    if gene_folders_with_input:
        ready_genes = [g for g in gene_folders_with_input if g['ready']]
        not_ready_genes = [g for g in gene_folders_with_input if not g['ready']]
        
        if ready_genes:
            st.success(f"✅ Found {len(ready_genes)} gene(s) **ready** for batch analysis")
            
            # Display ready genes
            for gene_info in ready_genes:
                st.info(f"**{gene_info['gene']}**: {gene_info['edited_count']} edited sample(s) ✅")
        
        if not_ready_genes:
            st.warning(f"⚠️ Found {len(not_ready_genes)} gene(s) with **missing files**")
            
            for gene_info in not_ready_genes:
                missing = ", ".join(gene_info['missing'])
                st.markdown(f"🔴 **{gene_info['gene']}**: {gene_info['edited_count']} edited samples, missing: `{missing}`")
        
        # Options
        col1, col2, col3 = st.columns(3)
        
        with col1:
            skip_genes = st.multiselect("Skip genes (optional)", 
                                       [g['gene'] for g in gene_folders_with_input],
                                       default=[g['gene'] for g in not_ready_genes],  # Auto-skip genes with issues
                                       key="batch_skip_genes")
        
        with col2:
            archive_previous = st.checkbox("Archive previous results", value=True, key="batch_archive")
            
        with col3:
            r_squared_correction = st.checkbox("Apply R² Correction", value=True, key="batch_r_squared")
        
        # Run batch analysis
        if ready_genes:
            if st.button("🚀 Run Batch Analysis", type="primary", width="stretch", key="batch_run_button"):
                run_batch_analysis(skip_genes, archive_previous, r_squared_correction)
        else:
            st.error("❌ No genes are ready for batch analysis. Please add missing files.")
        
        # Always display results if they exist
        display_batch_results()
    else:
        st.warning("No genes found with input files ready for batch analysis")
        
        # Show all gene folders and what's missing
        all_gene_folders = []
        for gene_dir in os.listdir(data_dir):
            gene_path = os.path.join(data_dir, gene_dir)
            if os.path.isdir(gene_path) and not gene_dir.startswith('.'):
                all_gene_folders.append(gene_dir)
        
        if all_gene_folders:
            st.markdown("### 📂 Gene Folders Found (but not ready)")
            
            for gene_dir in all_gene_folders:
                gene_path = os.path.join(data_dir, gene_dir)
                input_path = os.path.join(gene_path, "input")
                
                issues = []
                if not os.path.exists(input_path):
                    issues.append("No input/ folder")
                else:
                    if not os.path.exists(os.path.join(input_path, "control.ab1")):
                        issues.append("No control.ab1")
                    edited_count = len([f for f in os.listdir(input_path) 
                                       if f.startswith("edited") and f.endswith(".ab1")])
                    if edited_count == 0:
                        issues.append("No edited*.ab1 files")
                
                if not os.path.exists(os.path.join(gene_path, "grna.txt")):
                    issues.append("No grna.txt")
                if not os.path.exists(os.path.join(gene_path, "mrna.txt")):
                    issues.append("No mrna.txt")
                
                if issues:
                    st.markdown(f"🔴 **{gene_dir.upper()}**: {', '.join(issues)}")
                else:
                    st.markdown(f"✅ **{gene_dir.upper()}**: Ready")
        
        st.markdown("---")
        st.info("""
        **To prepare for batch analysis, each gene folder needs:**
        ```
        <data_dir>/<gene>/
        ├── grna.txt          # gRNA sequences (one per line)
        ├── mrna.txt          # mRNA/cDNA sequence
        └── input/
            ├── control.ab1   # Control sample AB1 file
            ├── editedA1.ab1  # Edited sample 1
            ├── editedA2.ab1  # Edited sample 2
            └──...           # More edited samples (edited*.ab1)
        ```
        """)


def check_gene_folder_requirements(gene_folder, gene_name):
    """Check what's missing in a gene folder for batch analysis."""
    issues = []
    
    logger.info(f"Checking requirements for {gene_name} at {gene_folder}")
    
    input_folder = os.path.join(gene_folder, "input")
    if not os.path.exists(input_folder):
        issues.append("No 'input/' folder")
        logger.warning(f"  No input folder: {input_folder}")
        return issues, None, []
    
    # Check for control file
    control_file = os.path.join(input_folder, "control.ab1")
    if not os.path.exists(control_file):
        issues.append("Missing control.ab1")
        logger.warning(f"  No control file: {control_file}")
        control_file = None
    else:
        logger.info(f"  Found control: {control_file}")
    
    # Check for edited files
    import glob
    edited_pattern = os.path.join(input_folder, "edited*.ab1")
    edited_files = glob.glob(edited_pattern)
    logger.info(f"  Looking for pattern: {edited_pattern}")
    logger.info(f"  Found {len(edited_files)} edited files")
    
    if not edited_files:
        issues.append("No edited*.ab1 files")
        # List what IS in the folder for debugging
        try:
            all_files = os.listdir(input_folder)
            logger.info(f"  Files in input folder: {all_files}")
        except Exception as e:
            logger.error(f"  Cannot list input folder: {e}")
    
    # Check for grna.txt
    grna_file = os.path.join(gene_folder, "grna.txt")
    if not os.path.exists(grna_file):
        issues.append("Missing grna.txt")
        logger.warning(f"  No grna file: {grna_file}")
    else:
        logger.info(f"  Found grna.txt")
    
    # Check for mrna.txt
    mrna_file = os.path.join(gene_folder, "mrna.txt")
    if not os.path.exists(mrna_file):
        issues.append("Missing mrna.txt")
        logger.warning(f"  No mrna file: {mrna_file}")
    else:
        logger.info(f"  Found mrna.txt")
    
    logger.info(f"  Issues found: {issues if issues else 'None'}")
    
    return issues, control_file, edited_files


def run_batch_analysis(skip_genes, archive_previous, r_squared_correction=True):
    """Run batch analysis on all prepared genes and store in session state."""
    from src.indel_analysis_multi import process_gene_folder
    
    # Create a simple args object
    class Args:
        def __init__(self):
            self.skip_genes = [g.lower() for g in skip_genes]
            self.force = False
            self.no_archive = not archive_previous
            self.r_squared_correction = r_squared_correction
    
    args = Args()
    
    # Process each gene
    data_dir = st.session_state.get('data_directory', '')
    if not data_dir or not os.path.isdir(data_dir):
        st.error("⚠️ Please set a valid data directory in the sidebar settings")
        return
    
    gene_folders = []
    
    for gene in os.listdir(data_dir):
        if gene.lower() not in args.skip_genes:
            gene_path = os.path.join(data_dir, gene)
            if os.path.isdir(gene_path) and os.path.exists(os.path.join(gene_path, "input")):
                gene_folders.append((gene_path, gene))
    
    if not gene_folders:
        st.warning("No gene folders found for batch analysis")
        return
    
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    # Store all results for display
    all_results = []
    successful = 0
    
    for i, (folder_path, gene_name) in enumerate(gene_folders):
        status_text.text(f"Processing {gene_name.upper()}...")
        
        # First check requirements
        issues, control_file, edited_files = check_gene_folder_requirements(folder_path, gene_name)
        
        if issues:
            all_results.append({
                'gene': gene_name.upper(),
                'folder': folder_path,
                'success': False,
                'error': "; ".join(issues)
            })
            progress_bar.progress((i + 1) / len(gene_folders))
            continue
        
        try:
            logger.info(f"Starting process_gene_folder for {gene_name}")
            result = process_gene_folder(folder_path, gene_name, args)
            logger.info(f"process_gene_folder returned: {result}")
            
            if result:
                successful += 1
                all_results.append({
                    'gene': gene_name.upper(),
                    'folder': folder_path,
                    'success': True,
                    'num_samples': len(edited_files)
                })
            else:
                all_results.append({
                    'gene': gene_name.upper(),
                    'folder': folder_path,
                    'success': False,
                    'error': "Processing returned False - check terminal for details"
                })
        except Exception as e:
            import traceback
            error_trace = traceback.format_exc()
            all_results.append({
                'gene': gene_name.upper(),
                'folder': folder_path,
                'success': False,
                'error': str(e)
            })
            logger.error(f"Batch analysis error for {gene_name}: {error_trace}")
            print(f"\n[ERROR] Batch analysis failed for {gene_name}:")
            print(error_trace)
        
        progress_bar.progress((i + 1) / len(gene_folders))
    
    status_text.text("Batch analysis complete!")
    st.success(f"✅ Successfully processed {successful}/{len(gene_folders)} genes")
    
    # Store results in session state for persistent display
    st.session_state.batch_results = all_results
    st.session_state.batch_analysis_done = True


def display_batch_results():
    """Display batch analysis results from session state."""
    if not st.session_state.get('batch_analysis_done', False):
        return
    
    all_results = st.session_state.get('batch_results', [])
    if not all_results:
        return
    
    # Display results with pagination
    st.markdown("---")
    st.markdown("## 📊 Batch Analysis Results")
    
    # Get successful results
    successful_results = [r for r in all_results if r['success']]
    
    if successful_results:
        # Select gene to view
        gene_options = [r['gene'] for r in successful_results]
        selected_gene = st.selectbox("Select Gene to View", gene_options, key="batch_gene_select")
        
        # Find selected gene's data
        current_result = next((r for r in successful_results if r['gene'] == selected_gene), None)
        
        if current_result:
            st.markdown(f"### Gene: {current_result['gene']}")
            
            # Find and display output files
            output_dir = os.path.join(current_result['folder'], "output")
            if os.path.exists(output_dir):
                # Find JSON results file
                json_files = [f for f in os.listdir(output_dir) if f.endswith('.json')]
                png_files = sorted([f for f in os.listdir(output_dir) if f.endswith('.png')])
                
                # Load and display results if JSON exists
                if json_files:
                    latest_json = sorted(json_files)[-1]
                    json_path = os.path.join(output_dir, latest_json)
                    
                    try:
                        with open(json_path, 'r') as f:
                            results_data = json.load(f)
                        
                        # Display summary metrics
                        samples = results_data.get('samples', [])
                        if samples:
                            st.markdown("#### 📋 Sample Results")
                            
                            # If samples is a list, iterate through it
                            if isinstance(samples, list):
                                for sample_data in samples:
                                    sample_name = sample_data.get('sample_name', 'Unknown')
                                    with st.expander(f"🧬 {sample_name}", expanded=False):
                                        col1, col2, col3, col4 = st.columns(4)
                                        with col1:
                                            eff = sample_data.get('editing_efficiency', 0)
                                            st.metric("Efficiency", f"{eff:.1f}%")
                                        with col2:
                                            wt = sample_data.get('wt_fraction', 100 - eff)
                                            st.metric("Wild-type", f"{wt:.1f}%")
                                        with col3:
                                            qual = sample_data.get('quality_score', 0)
                                            st.metric("R² Quality", f"{qual:.1f}%")
                                        with col4:
                                            conf = sample_data.get('confidence', 'Unknown')
                                            st.metric("Confidence", str(conf)[:15])
                                        
                                        # Display indel spectrum with TIDE-style colors (simple version, no nested expanders)
                                        indel_spectrum = sample_data.get('indel_spectrum', {})
                                        if indel_spectrum:
                                            display_indel_spectrum_simple(st, sample_data)
                            
                            # If samples is a dict, iterate through items
                            elif isinstance(samples, dict):
                                for sample_name, sample_data in samples.items():
                                    with st.expander(f"🧬 {sample_name}", expanded=False):
                                        col1, col2, col3, col4 = st.columns(4)
                                        with col1:
                                            eff = sample_data.get('editing_efficiency', 0)
                                            st.metric("Efficiency", f"{eff:.1f}%")
                                        with col2:
                                            wt = sample_data.get('wt_fraction', 100 - eff)
                                            st.metric("Wild-type", f"{wt:.1f}%")
                                        with col3:
                                            qual = sample_data.get('quality_score', 0)
                                            st.metric("R² Quality", f"{qual:.1f}%")
                                        with col4:
                                            conf = sample_data.get('confidence', 'Unknown')
                                            st.metric("Confidence", str(conf)[:15])
                        
                    except Exception as e:
                        st.warning(f"Could not load results: {e}")
                
                # Display images
                if png_files:
                    st.markdown("#### 🖼️ Analysis Plots")
                    
                    # Group images by type
                    summary_imgs = [f for f in png_files if 'summary' in f.lower()]
                    sample_imgs = [f for f in png_files if 'summary' not in f.lower()]
                    
                    # Show summary first
                    if summary_imgs:
                        st.markdown("##### Summary Plot")
                        for img_file in summary_imgs:
                            img_path = os.path.join(output_dir, img_file)
                            st.image(img_path, caption="Summary")
                    
                    # Show individual samples with selectbox for pagination
                    if sample_imgs:
                        st.markdown("##### Individual Sample Plots")
                        
                        # Use a unique key based on gene to avoid conflicts
                        img_select_key = f"sample_img_{selected_gene}"
                        
                        # Get current index from session state or default to 0
                        if img_select_key not in st.session_state:
                            st.session_state[img_select_key] = 0
                        
                        selected_idx = st.selectbox(
                            "Select Sample Plot", 
                            range(len(sample_imgs)),
                            format_func=lambda x: sample_imgs[x],
                            key=img_select_key
                        )
                        
                        selected_sample_img = sample_imgs[selected_idx]
                        img_path = os.path.join(output_dir, selected_sample_img)
                        st.image(img_path, caption=selected_sample_img)
                        
                        # Download button with unique key
                        with open(img_path, 'rb') as f:
                            st.download_button(
                                label="📥 Download This Plot",
                                data=f.read(),
                                file_name=selected_sample_img,
                                mime="image/png",
                                key=f"dl_{selected_gene}_{selected_idx}"
                            )
                        
                        # Download all as zip
                        st.markdown("---")
                        if st.button("📦 Download All Plots (ZIP)", key=f"dl_zip_{selected_gene}"):
                            import zipfile
                            import io
                            
                            zip_buffer = io.BytesIO()
                            with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zf:
                                for img_file in png_files:
                                    img_path = os.path.join(output_dir, img_file)
                                    zf.write(img_path, img_file)
                            
                            st.download_button(
                                label="📥 Download ZIP",
                                data=zip_buffer.getvalue(),
                                file_name=f"indel_analysis_{selected_gene}.zip",
                                mime="application/zip",
                                key=f"dl_zip_btn_{selected_gene}"
                            )
            else:
                st.info(f"Output directory not found for {current_result['gene']}")
    
    # Summary table of all genes
    st.markdown("---")
    st.markdown("### 📋 All Genes Summary")
    
    summary_data = []
    for r in all_results:
        status = '✅ Success' if r['success'] else '❌ Failed'
        error_msg = r.get('error', '') if not r['success'] else ''
        samples = r.get('num_samples', '-') if r['success'] else '-'
        
        summary_data.append({
            'Gene': r['gene'],
            'Status': status,
            'Samples': samples,
            'Issue': error_msg[:50] + '...' if len(error_msg) > 50 else error_msg,
            'Folder': r['folder']
        })

    st.dataframe(pd.DataFrame(summary_data), hide_index=True, width="stretch")
    
    # Show detailed errors for failed genes
    failed_genes = [r for r in all_results if not r['success']]
    if failed_genes:
        st.markdown("---")
        st.markdown("### ❌ Failed Gene Details")
        
        for r in failed_genes:
            with st.expander(f"🔴 {r['gene']} - {r.get('error', 'Unknown error')[:40]}...", expanded=False):
                st.markdown(f"**Folder:** `{r['folder']}`")
                st.markdown(f"**Error:** {r.get('error', 'Unknown error')}")
                
                # Show what's needed
                st.markdown("**Required files for batch analysis:**")
                st.markdown("""
                ```
                gene_folder/
                ├── grna.txt          # gRNA sequences (one per line)
                ├── mrna.txt          # mRNA/cDNA sequence
                └── input/
                    ├── control.ab1   # Control sample AB1 file
                    ├── editedA1.ab1  # Edited sample 1
                    ├── editedA2.ab1  # Edited sample 2
                    └──...           # More edited samples
                ```
                """)
    
    # Button to clear results and run again
    if st.button("🔄 Clear Results & Run New Analysis", key="clear_batch"):
        st.session_state.batch_analysis_done = False
        st.session_state.batch_results = []
        st.rerun()


def documentation_tab():
    """Documentation and help tab."""
    st.header("📚 Documentation")
    
    # Add new features highlight
    st.markdown("""
    <div class="success-box">
    🎉 <strong>New Enhanced Features!</strong><br>
    • Advanced primer design with database integration<br>
    • Comprehensive quality metrics and analysis<br>
    • Multiple output formats (report, order-ready, JSON)<br>
    • Real-time NCBI, Ensembl, and UniProt database checking<br>
    • Professional-grade primer analysis tools
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown("""
    ### Quick Start Guide
    
    #### 1. 🔍 gRNA Design
    - Enter your target gene name or Ensembl ID
    - Select species and Cas type (SpCas9, SaCas9, Cas12a, etc.)
    - Adjust filters (GC content, homopolymers, etc.)
    - Click "Find gRNAs" to get optimized guide sequences
    - Select and save gRNAs for downstream analysis
    
    #### 2. 🧪 Enhanced Primer Design
    **New Advanced Features:**
    - **Database Integration**: Automatic search of NCBI, Ensembl, UniProt
    - **Comprehensive Analysis**: Hairpin formation, primer dimers, specificity
    - **Quality Metrics**: Detailed scoring and recommendations
    - **Multiple Formats**: Full report, ready-to-order, JSON data
    - **Species Support**: Human, mouse, rat, and many others
    
    **Usage:**
    - Select gene and species from comprehensive list
    - Paste sequence or use local files
    - Provide gRNA sequences for optimal primer placement
    - Get two primer sets: PCR I (genomic) and PCR II (sequencing)
    - Download results in your preferred format
    
    #### 3. 📊 Indel Analysis
    - **Single Sample**: Upload control and edited AB1 files
    - **Batch Analysis**: Process multiple samples automatically
    - View editing efficiency and sequence traces
    - Download results and plots
    - Confidence scoring for result reliability
    
    ### Enhanced Primer Design Details
    
    #### Database Integration
    - **NCBI**: Gene information, sequence validation
    - **Ensembl**: Gene structure, variants, conservation
    - **UniProt**: Protein domains and functions (planned)
    
    #### Advanced Analysis Features
    - **Thermodynamic Analysis**: Tm calculation with salt correction
    - **Secondary Structure**: Hairpin and dimer prediction
    - **Quality Scoring**: Comprehensive primer evaluation
    - **Variant Checking**: Known SNPs in primer regions
    - **GC Distribution**: 5' and 3' region analysis
    
    #### Output Formats
    1. **Full Report**: Comprehensive analysis with all details
    2. **Ready-to-Order**: Copy-paste format for primer companies
    3. **JSON Data**: Machine-readable for further processing
    
    ### File Organization
    ```
    <your_data_directory>/      # Set in Settings panel
    ├── gene_name/
    │   ├── grna.txt            # gRNA sequences
    │   ├── mrna.txt            # mRNA sequence
    │   ├── primer_recommendations.txt
    │   ├── input/              # AB1 files for analysis
    │   │   ├── control.ab1
    │   │   └── edited*.ab1
    │   └── output/             # Analysis results
    ```
    
    **Important**: Data and code are SEPARATED. Set your data directory
    in the sidebar Settings panel or via CRISPR_DATA_DIR environment variable.
    
    ### Best Practices
    
    #### For gRNA Design
    - Use official gene symbols or Ensembl IDs
    - Check multiple Cas systems if needed
    - Consider off-target effects
    - Verify gRNA presence in target sequence
    
    #### For Primer Design
    - **Always validate with NCBI Primer-BLAST** before ordering
    - Order 2-3 primer sets for redundancy
    - Check for known variants in primer regions
    - Use gradient PCR to optimize annealing temperature
    - Consider genomic vs. cDNA differences
    
    #### For Indel Analysis
    - Use high-quality Sanger sequencing data
    - Ensure proper peak spacing and signal strength
    - Check confidence scores in results
    - Re-sequence if fallback method is used
    
    ### Comparison with Other Tools
    
    | Feature | Our Tool | Primer3 | IDT PrimerQuest | NEB Design |
    |---------|----------|---------|-----------------|------------|
    | Database Integration | ✅ | ❌ | Limited | ❌ |
    | Variant Checking | ✅ | ❌ | ✅ | ❌ |
    | CRISPR-Specific | ✅ | ❌ | Limited | ❌ |
    | Quality Analysis | ✅ | Basic | Good | Basic |
    | Multiple Formats | ✅ | ❌ | ❌ | ❌ |
    | Batch Processing | ✅ | ✅ | Limited | ❌ |
    
    ### Command Line Usage
    The GUI preserves all command-line functionality with enhancements:
    ```bash
    # Enhanced gRNA design
    python find_grna.py GENE_NAME --species human --cas-type SpCas9
    
    # Enhanced primer design with database checking
    python design_primers.py GENE_NAME --species human --check-databases
    
    # Multiple output formats
    python design_primers.py GENE_NAME --output-format detailed
    
    # Run indel analysis
    python run_analysis.py
    ```
    
    ### Troubleshooting
    
    #### gRNA Design Issues
    - **No gRNAs found**: Check gene name spelling, try Ensembl ID
    - **Low scores**: Adjust GC content range or other filters
    - **Off-target concerns**: Use stricter specificity settings
    
    #### Primer Design Issues
    - **No primers designed**: Check sequence length (min 100bp)
    - **High penalty scores**: Adjust Tm or GC content ranges
    - **Database search fails**: Check internet connection, try local files
    
    #### Indel Analysis Issues
    - **Low confidence**: Check AB1 file quality, re-sequence if needed
    - **Fallback method used**: Signal quality issues, consider new samples
    - **No cut site detected**: Verify gRNA sequences and positions
    
    ### Support & Updates
    - **GitHub**: Report issues and request features
    - **Documentation**: See [our GitHub documentation](https://github.com/InnovationLine/CasPINS/tree/main/docs) for detailed guides
    - **Updates**: Tool is actively maintained with regular enhancements
    """)
    
    # Add links to detailed documentation
    with st.expander("📖 View Detailed Documentation Files"):
        doc_files = [
            ("docs/grna_design_system.md", "gRNA Design System"),
            ("docs/primer_design_module.md", "Enhanced Primer Design Module"),
            ("docs/workflow_guide.md", "Complete Workflow Guide"),
            ("docs/adding_new_genes.md", "Adding New Genes"),
            ("docs/fallback_method_guide.md", "Fallback Method Guide"),
            ("docs/modular_architecture.md", "System Architecture")
        ]
        
        for doc_file, title in doc_files:
            github_url = f"https://github.com/InnovationLine/CasPINS/blob/main/{doc_file}"
            st.markdown(f"#### 📄 [{title}]({github_url})")
            if os.path.exists(doc_file):
                with st.expander(f"Read {title} inline"):
                    with open(doc_file, 'r', encoding='utf-8') as f:
                        st.markdown(f.read())
            st.markdown("---")
    
    # Add version information
    st.markdown("---")
    st.markdown("""
    ### Version Information
    - **CasPINS - Cas-Primer-Indel Suite**: v2.0 (Enhanced)
    - **Primer Design**: v2.0 with database integration
    - **gRNA Design**: v1.5 with multi-Cas support
    - **Indel Analysis**: v1.3 with confidence scoring
    
    **Latest Updates:**
    - ✨ Enhanced primer design with professional-grade analysis
    - 🔗 Direct database integration (NCBI, Ensembl, UniProt)
    - 📊 Comprehensive quality metrics and recommendations
    - 📁 Multiple output formats for different use cases
    - 🧬 Expanded species support and Cas system options
    """)
    
    # Add contact information
    with st.expander("📞 Contact & Support"):
        st.markdown("""
        **For Questions or Issues:**
        - GitHub Issues: Report bugs or request features
        - Email: crispr_analysis@example.com
        - Documentation: See `/docs/` folder for detailed guides
        
        **Contributing:**
        - Fork the repository on GitHub
        - Submit pull requests for improvements
        - Report bugs and suggest features
        
        **Citation:**
        If you use this tool in your research, please cite:
        "CasPINS - Cas-Primer-Indel Suite: An integrated platform for gRNA design, primer design, and indel analysis"
        """)


def cli_entry():
    """Command-line entry point."""
    print("Starting CasPINS - Cas-Primer-Indel Suite GUI...")
    print("The GUI will open in your default web browser.")
    print("Press Ctrl+C to stop the server.")
    
    # Run streamlit
    sys.argv = ["streamlit", "run", __file__]
    sys.exit(st.cli.main())


if __name__ == "__main__":
    main() 