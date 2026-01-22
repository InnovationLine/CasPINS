"""
Main gRNA Designer Module
Integrates all components for comprehensive gRNA design

Supports:
- Multiple input types: Gene symbols, Ensembl IDs, RefSeq IDs, Genomic coordinates
- Multiple editing modes: Knockout, Knock-in (HDR), CRISPRa, CRISPRi, Base editing
- Genomic location output with chr:start-end format
- Multiple CRISPR systems and species
"""

from typing import List, Dict, Optional, Union, Tuple
import json
import os
import re
import logging
from datetime import datetime
from enum import Enum

# Configure logger
logger = logging.getLogger(__name__)

from .core.grna_generator import GRNAGenerator
from .core.sequence_analyzer import SequenceAnalyzer
from .scoring.scoring_engine import ScoringEngine
from .database.genome_manager import GenomeManager


class EditingMode(Enum):
    """CRISPR editing modes."""
    KNOCKOUT = "knockout"              # Standard gene disruption via NHEJ
    KNOCKIN_HDR = "knockin_hdr"        # Precise insertion via HDR
    ACTIVATION = "crispra"             # CRISPRa - Gene activation (dCas9-VP64, etc.)
    REPRESSION = "crispri"             # CRISPRi - Gene repression (dCas9-KRAB, etc.)
    BASE_EDITING_CBE = "base_edit_cbe" # Cytosine base editor (C>T)
    BASE_EDITING_ABE = "base_edit_abe" # Adenine base editor (A>G)
    PRIME_EDITING = "prime_edit"       # Prime editing with pegRNA
    NANOPORE_ENRICHMENT = "nanopore"   # Cas9-guided nanopore enrichment


class TargetType(Enum):
    """Types of target input."""
    GENE_SYMBOL = "gene_symbol"        # TP53, DDC, BRCA1
    ENSEMBL_GENE = "ensembl_gene"      # ENSG00000107669
    ENSEMBL_TRANSCRIPT = "ensembl_tx"  # ENST00000358428
    REFSEQ_MRNA = "refseq_mrna"        # NM_003054.5
    REFSEQ_GENE = "refseq_gene"        # NG_008849.1
    GENOMIC_COORDS = "genomic_coords"  # chr8:20100000-20150000
    SEQUENCE = "sequence"              # Raw DNA sequence


class GRNADesigner:
    """
    Main class for designing guide RNAs.
    
    Supports multiple input formats and editing modes.
    """
    
    # Editing mode configurations
    EDITING_MODE_CONFIG = {
        EditingMode.KNOCKOUT: {
            'target_region': 'exons',
            'description': 'Gene disruption via NHEJ',
            'grna_position': 'early_exons',
            'notes': 'Target early constitutive exons for complete knockout'
        },
        EditingMode.KNOCKIN_HDR: {
            'target_region': 'specific',
            'description': 'Precise insertion via HDR',
            'grna_position': 'near_target',
            'notes': 'Cut site should be <10bp from desired insertion point'
        },
        EditingMode.ACTIVATION: {
            'target_region': 'promoter',
            'description': 'CRISPRa - Gene activation',
            'grna_position': '200bp_upstream_tss',
            'notes': 'Target -200 to -50 upstream of TSS for optimal activation'
        },
        EditingMode.REPRESSION: {
            'target_region': 'promoter',
            'description': 'CRISPRi - Gene repression',
            'grna_position': 'around_tss',
            'notes': 'Target +50 to -50 around TSS for optimal repression'
        },
        EditingMode.BASE_EDITING_CBE: {
            'target_region': 'specific',
            'description': 'C>T (or G>A) base editing',
            'grna_position': 'editing_window',
            'notes': 'Target C must be in positions 4-8 of protospacer (activity window)'
        },
        EditingMode.BASE_EDITING_ABE: {
            'target_region': 'specific',
            'description': 'A>G (or T>C) base editing',
            'grna_position': 'editing_window',
            'notes': 'Target A must be in positions 4-7 of protospacer (activity window)'
        },
        EditingMode.PRIME_EDITING: {
            'target_region': 'specific',
            'description': 'Precise editing with pegRNA',
            'grna_position': 'near_target',
            'notes': 'Nick site should be 3bp from desired edit'
        },
        EditingMode.NANOPORE_ENRICHMENT: {
            'target_region': 'flanking',
            'description': 'Cas9-guided nanopore enrichment',
            'grna_position': 'flanking_region',
            'notes': 'Design gRNAs flanking region of interest for enrichment'
        }
    }
    
    def __init__(self, species: str = 'human', assembly: str = None,
                 cas_type: str = 'SpCas9', cache_dir: str = 'cache',
                 editing_mode: EditingMode = EditingMode.KNOCKOUT):
        """
        Initialize gRNA designer.
        
        Args:
            species: Target species (e.g., 'human', 'mouse', 'homo_sapiens')
            assembly: Genome assembly version (e.g., 'GRCh38', 'hg38')
            cas_type: CRISPR-Cas system type (e.g., 'SpCas9', 'SaCas9', 'Cas12a')
            cache_dir: Directory for caching data
            editing_mode: Type of editing (knockout, activation, repression, etc.)
        """
        self.species = self._normalize_species(species)
        self.assembly = assembly
        self.cas_type = cas_type
        self.editing_mode = editing_mode
        
        logger.info(f"GRNADesigner initialized: species={self.species}, cas={cas_type}")
        
        # Initialize components
        self.generator = GRNAGenerator(cas_type)
        self.analyzer = SequenceAnalyzer(os.path.join(cache_dir, 'sequences'))
        self.scorer = ScoringEngine()
        self.genome_manager = GenomeManager(os.path.join(cache_dir, 'genomes'))
        
        # Design history
        self.history_file = os.path.join(cache_dir, 'design_history.json')
    
    def _normalize_species(self, species: str) -> str:
        """Normalize species name to Ensembl format (genus_species)."""
        # Comprehensive mapping from common names and display formats to Ensembl names
        species_map = {
            # Common names
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
            'c_elegans': 'caenorhabditis_elegans',
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
            'ferret': 'mustela_putorius_furo',
            'mustela putorius furo': 'mustela_putorius_furo',
        }
        
        # Clean up input - handle display format "Human (Homo sapiens)"
        species_clean = species.split('(')[0].strip().lower().replace('_', ' ')
        
        # Also check parenthetical content
        paren_content = ''
        if '(' in species and ')' in species:
            paren_content = species.split('(')[1].split(')')[0].strip().lower()
        
        # Try exact match first
        if species_clean in species_map:
            return species_map[species_clean]
        
        if paren_content and paren_content in species_map:
            return species_map[paren_content]
        
        # Try with underscores
        species_underscore = species_clean.replace(' ', '_')
        if species_underscore in species_map:
            return species_map[species_underscore]
        
        # Try partial match
        for key, value in species_map.items():
            if key in species_clean or species_clean in key:
                return value
        
        # Default: convert to genus_species format
        return species_clean.replace(' ', '_')
    
    def _detect_target_type(self, target: str) -> TargetType:
        """
        Auto-detect the type of target input.
        
        Args:
            target: Input string to analyze
            
        Returns:
            TargetType enum value
        """
        target = target.strip()
        
        # Genomic coordinates: chr1:12345-67890 or chr1:12345..67890
        coord_pattern = r'^chr[\dXYMT]+[:\s]+\d+[-\.]{1,2}\d+$'
        if re.match(coord_pattern, target, re.IGNORECASE):
            return TargetType.GENOMIC_COORDS
        
        # Ensembl gene ID: ENSG00000123456
        if re.match(r'^ENS[A-Z]*G\d{11}$', target, re.IGNORECASE):
            return TargetType.ENSEMBL_GENE
        
        # Ensembl transcript ID: ENST00000123456
        if re.match(r'^ENS[A-Z]*T\d{11}$', target, re.IGNORECASE):
            return TargetType.ENSEMBL_TRANSCRIPT
        
        # RefSeq mRNA: NM_123456.1 or NM_123456
        if re.match(r'^NM_\d+(\.\d+)?$', target, re.IGNORECASE):
            return TargetType.REFSEQ_MRNA
        
        # RefSeq gene: NG_123456.1
        if re.match(r'^NG_\d+(\.\d+)?$', target, re.IGNORECASE):
            return TargetType.REFSEQ_GENE
        
        # Raw sequence (only ATCGN characters, >20bp)
        if re.match(r'^[ATCGNatcgn]{20,}$', target):
            return TargetType.SEQUENCE
        
        # Default: treat as gene symbol
        return TargetType.GENE_SYMBOL
    
    def _parse_genomic_coordinates(self, coord_string: str) -> Tuple[str, int, int]:
        """
        Parse genomic coordinates from string.
        
        Args:
            coord_string: String like 'chr8:20100000-20150000'
            
        Returns:
            Tuple of (chromosome, start, end)
        """
        # Normalize format
        coord_string = coord_string.strip().replace('..', '-').replace(' ', '')
        
        # Parse
        match = re.match(r'^(chr[\dXYMT]+)[:\s]+(\d+)-(\d+)$', coord_string, re.IGNORECASE)
        if match:
            chrom = match.group(1)
            start = int(match.group(2))
            end = int(match.group(3))
            return (chrom, start, end)
        
        raise ValueError(f"Invalid coordinate format: {coord_string}")
    
    def _format_genomic_location(self, chromosome: str, start: int, end: int, 
                                  strand: int = 1, assembly: str = None) -> str:
        """
        Format genomic location in standard notation.
        
        Args:
            chromosome: Chromosome name
            start: Start position
            end: End position
            strand: Strand (1 or -1)
            assembly: Assembly version
            
        Returns:
            Formatted string like 'chr8:20,100,000-20,150,000 (+) [GRCh38]'
        """
        # Format numbers with commas
        start_fmt = f"{start:,}"
        end_fmt = f"{end:,}"
        
        # Strand symbol
        strand_sym = "(+)" if strand == 1 else "(-)"
        
        # Construct string
        location = f"{chromosome}:{start_fmt}-{end_fmt} {strand_sym}"
        
        if assembly:
            location += f" [{assembly}]"
        
        return location
        
    def design_grnas(self, target: Union[str, Dict], 
                    n_results: int = 10,
                    target_type: str = 'auto',
                    filters: Dict = None,
                    include_off_targets: bool = True,
                    editing_mode: EditingMode = None) -> Dict:
        """
        Design gRNAs for a target.
        
        Supports multiple input types:
        - Gene symbols: TP53, DDC, BRCA1
        - Ensembl IDs: ENSG00000107669, ENST00000358428
        - RefSeq IDs: NM_003054.5, NG_008849.1
        - Genomic coordinates: chr8:20100000-20150000
        - Raw sequences: ATCGATCG...
        
        Args:
            target: Gene name, Ensembl ID, RefSeq ID, coordinates, or sequence
            n_results: Number of top gRNAs to return
            target_type: Type of target ('auto', 'gene', 'sequence', 'coordinates')
            filters: Filtering parameters
            include_off_targets: Whether to analyze off-targets
            editing_mode: Override editing mode for this design
            
        Returns:
            Dictionary with design results including genomic locations
        """
        # Use instance editing mode if not overridden
        if editing_mode is None:
            editing_mode = self.editing_mode
        
        # Auto-detect target type if needed
        if target_type == 'auto' and isinstance(target, str):
            detected_type = self._detect_target_type(target)
            target_type = detected_type.value
        
        results = {
            'target': target,
            'target_type': target_type,
            'cas_type': self.cas_type,
            'species': self.species,
            'assembly': self.assembly,
            'editing_mode': editing_mode.value if isinstance(editing_mode, EditingMode) else editing_mode,
            'timestamp': datetime.now().isoformat(),
            'grnas': [],
            'metadata': {}
        }
        
        # Get target sequence based on input type
        sequence_data = None
        
        if target_type in ['gene', 'gene_symbol', TargetType.GENE_SYMBOL.value]:
            sequence_data = self._get_gene_sequence(target)
        elif target_type in ['ensembl_gene', TargetType.ENSEMBL_GENE.value]:
            sequence_data = self._get_ensembl_sequence(target, 'gene')
        elif target_type in ['ensembl_tx', TargetType.ENSEMBL_TRANSCRIPT.value]:
            sequence_data = self._get_ensembl_sequence(target, 'transcript')
        elif target_type in ['refseq_mrna', TargetType.REFSEQ_MRNA.value]:
            sequence_data = self._get_refseq_sequence(target)
        elif target_type in ['genomic_coords', TargetType.GENOMIC_COORDS.value]:
            sequence_data = self._get_coordinate_sequence(target)
        elif target_type in ['sequence', TargetType.SEQUENCE.value]:
            sequence_data = target if isinstance(target, dict) else {'sequence': target}
        else:
            # Fallback: try as gene symbol
            sequence_data = self._get_gene_sequence(target)
        
        if not sequence_data:
            results['error'] = f"Could not retrieve sequence for {target}"
            return results
            
        # Store enhanced metadata
        results['metadata'] = {
            'gene_name': sequence_data.get('gene_name', target),
            'gene_symbol': sequence_data.get('gene_symbol', ''),
            'ensembl_id': sequence_data.get('ensembl_id', ''),
            'refseq_id': sequence_data.get('refseq_id', ''),
            'chromosome': sequence_data.get('chromosome'),
            'start': sequence_data.get('start'),
            'end': sequence_data.get('end'),
            'strand': sequence_data.get('strand'),
            'assembly': sequence_data.get('assembly', self.assembly),
            'source': sequence_data.get('source', 'unknown'),
            'sequence_length': len(sequence_data.get('sequence', ''))
        }
        
        # Add formatted genomic location
        if results['metadata']['chromosome'] and results['metadata']['start']:
            results['metadata']['genomic_location'] = self._format_genomic_location(
                results['metadata']['chromosome'],
                results['metadata']['start'],
                results['metadata']['end'],
                results['metadata']['strand'] or 1,
                results['metadata']['assembly']
            )
        
        # Add editing mode information
        if isinstance(editing_mode, EditingMode) and editing_mode in self.EDITING_MODE_CONFIG:
            results['editing_mode_info'] = self.EDITING_MODE_CONFIG[editing_mode]
        
        # Find all possible gRNAs
        all_grnas = self.generator.find_all_grnas(sequence_data['sequence'])
        print(f"Found {len(all_grnas)} potential gRNAs")
        
        # Apply editing mode-specific filtering
        if isinstance(editing_mode, EditingMode):
            all_grnas = self._filter_for_editing_mode(all_grnas, sequence_data, editing_mode)
            print(f"After editing mode filtering: {len(all_grnas)} gRNAs")
        
        # Apply user filters
        if filters:
            filtered_grnas = self.generator.comprehensive_filter(all_grnas, filters)
            print(f"After user filtering: {len(filtered_grnas)} gRNAs")
        else:
            filtered_grnas = all_grnas
            
        # Score gRNAs
        for grna in filtered_grnas:
            # Get sequence context if available
            if 'position' in grna:
                context = self.analyzer.get_sequence_context(
                    sequence_data['sequence'],
                    grna['position'],
                    upstream=30,
                    downstream=30
                )
                grna['context'] = context
                
            # Calculate scores
            grna['scores'] = self.scorer.calculate_all_scores(
                grna['sequence'],
                grna.get('context'),
                grna.get('pam', 'NGG')
            )
            
        # Rank by composite score
        ranked_grnas = self.scorer.rank_grnas(filtered_grnas)
        
        # Get top N results
        top_grnas = ranked_grnas[:n_results]
        
        # Add genomic coordinates and location strings for each gRNA
        sequence_start = sequence_data.get('start') or sequence_data.get('sequence_start') or 0
        chromosome = results['metadata'].get('chromosome', '')
        assembly = results['metadata'].get('assembly', '')
        strand = results['metadata'].get('strand', 1)
        
        for grna in top_grnas:
            if 'position' in grna:
                grna_start = sequence_start + grna['position']
                grna_end = grna_start + len(grna['sequence'])
                
                grna['genomic_start'] = grna_start
                grna['genomic_end'] = grna_end
                
                if chromosome:
                    grna['genomic_location'] = self._format_genomic_location(
                        chromosome, grna_start, grna_end, strand, assembly
                    )
                    # Also add simple format
                    grna['location_simple'] = f"{chromosome}:{grna_start}-{grna_end}"
        
        # Analyze off-targets if requested
        if include_off_targets and self.genome_manager.get_genome(self.species, self.assembly):
            for grna in top_grnas:
                grna['off_targets'] = self._analyze_off_targets(grna['sequence'])
                    
        results['grnas'] = top_grnas
        
        # Save to history
        self._save_to_history(results)
        
        return results
    
    def _filter_for_editing_mode(self, grnas: List[Dict], sequence_data: Dict, 
                                  editing_mode: EditingMode) -> List[Dict]:
        """
        Filter gRNAs based on editing mode requirements.
        
        Args:
            grnas: List of gRNA dictionaries
            sequence_data: Sequence data with metadata
            editing_mode: Editing mode to filter for
            
        Returns:
            Filtered list of gRNAs
        """
        if editing_mode == EditingMode.ACTIVATION:
            # CRISPRa: Target promoter region, -200 to -50 from TSS
            # For simplicity, target first 500bp of sequence
            return [g for g in grnas if g.get('position', 0) < 500]
        
        elif editing_mode == EditingMode.REPRESSION:
            # CRISPRi: Target around TSS, -50 to +300
            return [g for g in grnas if g.get('position', 0) < 800]
        
        elif editing_mode == EditingMode.BASE_EDITING_CBE:
            # CBE: Need C in positions 4-8 of protospacer
            filtered = []
            for g in grnas:
                seq = g['sequence']
                if len(seq) >= 8 and 'C' in seq[3:8]:
                    g['base_editing_targets'] = [i+1 for i, b in enumerate(seq[3:8]) if b == 'C']
                    filtered.append(g)
            return filtered
        
        elif editing_mode == EditingMode.BASE_EDITING_ABE:
            # ABE: Need A in positions 4-7 of protospacer
            filtered = []
            for g in grnas:
                seq = g['sequence']
                if len(seq) >= 7 and 'A' in seq[3:7]:
                    g['base_editing_targets'] = [i+1 for i, b in enumerate(seq[3:7]) if b == 'A']
                    filtered.append(g)
            return filtered
        
        # For other modes, return all gRNAs
        return grnas
    
    def _get_ensembl_sequence(self, ensembl_id: str, id_type: str = 'gene') -> Optional[Dict]:
        """
        Get sequence data from Ensembl ID.
        
        Args:
            ensembl_id: Ensembl gene or transcript ID
            id_type: 'gene' or 'transcript'
            
        Returns:
            Sequence data dictionary
        """
        # Use sequence analyzer to fetch from Ensembl
        try:
            sequence_data = self.analyzer.get_gene_sequence(
                ensembl_id,
                self.species,
                expand_5p=500,
                expand_3p=500
            )
            
            if sequence_data:
                sequence_data['source'] = 'ensembl'
                sequence_data['ensembl_id'] = ensembl_id
                
            return sequence_data
            
        except Exception as e:
            print(f"Error fetching Ensembl sequence: {e}")
            return None
    
    def _get_refseq_sequence(self, refseq_id: str) -> Optional[Dict]:
        """
        Get sequence data from RefSeq ID.
        
        Args:
            refseq_id: RefSeq accession (NM_*, NG_*, etc.)
            
        Returns:
            Sequence data dictionary
        """
        import requests
        
        try:
            # Use NCBI E-utilities to fetch sequence
            base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
            
            # First, get the sequence
            efetch_url = f"{base_url}/efetch.fcgi"
            params = {
                'db': 'nuccore',
                'id': refseq_id,
                'rettype': 'fasta',
                'retmode': 'text'
            }
            
            response = requests.get(efetch_url, params=params, timeout=30)
            
            if response.status_code == 200:
                fasta_text = response.text
                
                # Parse FASTA
                lines = fasta_text.strip().split('\n')
                header = lines[0] if lines else ''
                sequence = ''.join(lines[1:]).upper()
                
                # Extract gene name from header
                gene_name = refseq_id
                if '|' in header:
                    parts = header.split('|')
                    if len(parts) > 1:
                        gene_name = parts[-1].split()[0] if parts[-1] else refseq_id
                
                return {
                    'sequence': sequence,
                    'gene_name': gene_name,
                    'refseq_id': refseq_id,
                    'source': 'ncbi_refseq'
                }
                
        except Exception as e:
            print(f"Error fetching RefSeq sequence: {e}")
        
        return None
    
    def _get_coordinate_sequence(self, coord_string: str) -> Optional[Dict]:
        """
        Get sequence from genomic coordinates.
        
        Args:
            coord_string: Coordinates like 'chr8:20100000-20150000'
            
        Returns:
            Sequence data dictionary
        """
        try:
            chromosome, start, end = self._parse_genomic_coordinates(coord_string)
            
            # Fetch sequence from Ensembl
            import requests
            
            # Convert species name for Ensembl REST API
            species_map = {
                'homo_sapiens': 'human',
                'mus_musculus': 'mouse',
                'rattus_norvegicus': 'rat',
                'danio_rerio': 'zebrafish',
            }
            species_name = species_map.get(self.species, self.species)
            
            # Ensembl REST API for sequence by region
            url = f"https://rest.ensembl.org/sequence/region/{species_name}/{chromosome}:{start}..{end}"
            
            response = requests.get(url, headers={"Content-Type": "text/plain"}, timeout=30)
            
            if response.status_code == 200:
                sequence = response.text.strip().upper()
                
                return {
                    'sequence': sequence,
                    'gene_name': f"{chromosome}:{start}-{end}",
                    'chromosome': chromosome,
                    'start': start,
                    'end': end,
                    'strand': 1,
                    'assembly': self.assembly or 'GRCh38',
                    'source': 'ensembl_coordinates'
                }
                
        except Exception as e:
            print(f"Error fetching coordinate sequence: {e}")
        
        return None
    
    def design_batch(self, targets: List[str], 
                    target_type: str = 'gene',
                    **kwargs) -> List[Dict]:
        """
        Design gRNAs for multiple targets.
        
        Args:
            targets: List of target identifiers
            target_type: Type of targets
            **kwargs: Additional parameters for design_grnas
            
        Returns:
            List of design results
        """
        results = []
        
        for i, target in enumerate(targets):
            print(f"\nProcessing target {i+1}/{len(targets)}: {target}")
            result = self.design_grnas(target, target_type=target_type, **kwargs)
            results.append(result)
            
        return results
    
    def _get_gene_sequence(self, gene_id: str, data_dir: str = None) -> Optional[Dict]:
        """Get gene sequence data from Ensembl or local files."""
        # Try to get from sequence analyzer (Ensembl API)
        sequence_data = self.analyzer.get_gene_sequence(
            gene_id, 
            self.species,
            expand_5p=500,
            expand_3p=500
        )
        
        if not sequence_data and data_dir:
            # Try local files if data_dir is provided
            local_file = os.path.join(data_dir, gene_id.lower(), "mrna.txt")
            if os.path.exists(local_file):
                with open(local_file, 'r') as f:
                    sequence = f.read().strip()
                    
                sequence_data = {
                    'gene_name': gene_id,
                    'sequence': sequence,
                    'source': 'local_file'
                }
                
                # Load genomic metadata if available
                metadata_file = os.path.join(data_dir, "gene_metadata.json")
                if os.path.exists(metadata_file):
                    with open(metadata_file, 'r') as f:
                        metadata = json.load(f)
                    
                    if self.species in metadata and gene_id.upper() in metadata[self.species]:
                        gene_meta = metadata[self.species][gene_id.upper()]
                        sequence_data.update({
                            'chromosome': gene_meta.get('chromosome'),
                            'start': gene_meta.get('start'),
                            'end': gene_meta.get('end'),
                            'strand': gene_meta.get('strand'),
                            'ensembl_id': gene_meta.get('ensembl_id'),
                            'gene_symbol': gene_meta.get('gene_symbol'),
                            'description': gene_meta.get('description'),
                            'assembly': gene_meta.get('assembly', self.assembly)
                        })
                else:
                    # Try to fetch genomic data from Ensembl for any gene
                    print(f"Fetching genomic metadata for {gene_id} from Ensembl...")
                    genomic_data = self.analyzer._search_gene_by_symbol(gene_id, self.species)
                    if genomic_data:
                        sequence_data.update({
                            'chromosome': genomic_data.get('chromosome'),
                            'start': genomic_data.get('start'),
                            'end': genomic_data.get('end'),
                            'strand': genomic_data.get('strand'),
                            'ensembl_id': genomic_data.get('gene_id'),
                            'gene_symbol': gene_id,
                            'assembly': genomic_data.get('assembly_name')
                        })
                
        return sequence_data
    
    def _analyze_off_targets(self, grna_sequence: str, 
                           max_mismatches: int = 3) -> List[Dict]:
        """Analyze potential off-targets."""
        off_targets = []
        
        # Search genome for similar sequences
        matches = self.genome_manager.search_sequence(
            grna_sequence,
            self.species,
            max_mismatches=max_mismatches,
            assembly=self.assembly
        )
        
        # Calculate CFD scores
        for match in matches:
            if match['sequence'] != grna_sequence:  # Skip on-target
                cfd_score = self.scorer.calculate_cfd_score(
                    grna_sequence,
                    match['sequence']
                )
                
                match['cfd_score'] = cfd_score
                off_targets.append(match)
                
        # Sort by CFD score
        off_targets.sort(key=lambda x: x['cfd_score'], reverse=True)
        
        return off_targets[:10]  # Return top 10 off-targets
    
    def _save_to_history(self, results: Dict):
        """Save design results to history."""
        try:
            # Load existing history
            if os.path.exists(self.history_file):
                with open(self.history_file, 'r') as f:
                    history = json.load(f)
            else:
                history = []
                
            # Add new results
            history.append(results)
            
            # Save updated history
            os.makedirs(os.path.dirname(self.history_file), exist_ok=True)
            with open(self.history_file, 'w') as f:
                json.dump(history, f, indent=2)
                
        except Exception as e:
            print(f"Warning: Could not save to history: {str(e)}")
    
    def export_results(self, results: Dict, format: str = 'csv',
                      output_file: str = None) -> str:
        """
        Export results in various formats.
        
        Args:
            results: Design results dictionary
            format: Export format ('csv', 'json', 'genbank')
            output_file: Output filename
            
        Returns:
            Path to exported file
        """
        if not output_file:
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            gene_name = results['metadata'].get('gene_name', 'unknown')
            output_file = f"grna_design_{gene_name}_{timestamp}.{format}"
            
        if format == 'csv':
            self._export_csv(results, output_file)
        elif format == 'json':
            with open(output_file, 'w') as f:
                json.dump(results, f, indent=2)
        elif format == 'genbank':
            self._export_genbank(results, output_file)
        else:
            raise ValueError(f"Unsupported format: {format}")
            
        return output_file
    
    def _export_csv(self, results: Dict, output_file: str):
        """Export results as CSV."""
        import csv
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Write header
            writer.writerow([
                'Rank', 'gRNA Sequence', 'PAM', 'Position', 'Strand',
                'GC Content', 'Composite Score', 'Doench Score', 
                'Moreno-Mateos Score', 'Xu Score', 'Off-Target Count'
            ])
            
            # Write gRNA data
            for i, grna in enumerate(results['grnas']):
                writer.writerow([
                    i + 1,
                    grna['sequence'],
                    grna.get('pam', ''),
                    grna.get('position', ''),
                    grna.get('strand', ''),
                    grna.get('gc_content', ''),
                    grna['scores'].get('composite', ''),
                    grna['scores'].get('doench_2016', ''),
                    grna['scores'].get('moreno_mateos', ''),
                    grna['scores'].get('xu', ''),
                    len(grna.get('off_targets', []))
                ])
    
    def _export_genbank(self, results: Dict, output_file: str):
        """Export results as GenBank format."""
        from Bio import SeqRecord, SeqFeature
        from Bio.Seq import Seq
        
        # Create sequence record
        if 'sequence' in results.get('metadata', {}):
            sequence = results['metadata']['sequence']
        else:
            # Use a placeholder sequence
            sequence = "N" * 1000
            
        record = SeqRecord.SeqRecord(
            Seq(sequence),
            id=results['metadata'].get('gene_name', 'unknown'),
            description=f"gRNA design for {results['target']}"
        )
        
        # Add gRNA features
        for i, grna in enumerate(results['grnas']):
            if 'position' in grna:
                feature = SeqFeature.SeqFeature(
                    SeqFeature.FeatureLocation(
                        grna['position'],
                        grna['position'] + len(grna['sequence'])
                    ),
                    type="misc_feature",
                    qualifiers={
                        'label': f"gRNA_{i+1}",
                        'note': f"Score: {grna['scores']['composite']:.3f}",
                        'sequence': grna['sequence']
                    }
                )
                record.features.append(feature)
                
        # Write to file
        from Bio import SeqIO
        SeqIO.write(record, output_file, "genbank")
    
    def validate_grna(self, grna_sequence: str) -> Dict:
        """
        Validate a user-provided gRNA sequence.
        
        Args:
            grna_sequence: gRNA sequence to validate
            
        Returns:
            Validation results
        """
        validation = {
            'sequence': grna_sequence,
            'valid': True,
            'warnings': [],
            'errors': []
        }
        
        # Check length
        if len(grna_sequence) != 20:
            validation['warnings'].append(f"Non-standard length: {len(grna_sequence)} (expected 20)")
            
        # Check for valid bases
        valid_bases = set('ATCG')
        invalid_bases = set(grna_sequence.upper()) - valid_bases
        if invalid_bases:
            validation['errors'].append(f"Invalid bases: {invalid_bases}")
            validation['valid'] = False
            
        # Check GC content
        gc_content = (grna_sequence.count('G') + grna_sequence.count('C')) / len(grna_sequence) * 100
        if gc_content < 20 or gc_content > 80:
            validation['warnings'].append(f"Extreme GC content: {gc_content:.1f}%")
            
        # Check for homopolymers
        if self.generator.check_homopolymers(grna_sequence):
            validation['warnings'].append("Contains homopolymer runs")
            
        # Calculate scores
        if validation['valid']:
            validation['scores'] = self.scorer.calculate_all_scores(grna_sequence)
            
        return validation 