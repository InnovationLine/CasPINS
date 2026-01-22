"""
Sequence Analyzer Module
Handles gene sequence retrieval, analysis, and context extraction
"""

import requests
from typing import Dict, List, Tuple, Optional
from Bio import SeqIO
from Bio.Seq import Seq
import json
import os
import logging

# Configure logger
logger = logging.getLogger(__name__)


class SequenceAnalyzer:
    """Analyze gene sequences and extract relevant information."""
    
    # API endpoints for sequence retrieval
    ENSEMBL_REST = "https://rest.ensembl.org"
    NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    
    def __init__(self, cache_dir: str = "cache/sequences"):
        """
        Initialize sequence analyzer.
        
        Args:
            cache_dir: Directory for caching sequences
        """
        self.cache_dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)
        
    def get_gene_sequence(self, gene_id: str, species: str = "human",
                         expand_5p: int = 500, expand_3p: int = 500) -> Dict:
        """
        Retrieve gene sequence from Ensembl or NCBI.

        Args:
            gene_id: Gene identifier (Ensembl ID, gene symbol, or NCBI ID)
            species: Species name (e.g., 'homo_sapiens', 'rattus_norvegicus')
            expand_5p: Basepairs to include upstream
            expand_3p: Basepairs to include downstream

        Returns:
            Dictionary with sequence information
        """
        logger.info(f"Getting sequence for: {gene_id} (species: {species})")
        
        # Normalize species
        species_normalized = self._normalize_species_for_ensembl(species)
        logger.info(f"Normalized species: {species_normalized}")
        
        # Check cache first
        cache_file = os.path.join(self.cache_dir, f"{gene_id}_{species_normalized}.json")
        if os.path.exists(cache_file):
            logger.info(f"Found cached sequence: {cache_file}")
            with open(cache_file, 'r') as f:
                return json.load(f)

        # Try different retrieval methods
        result = None

        # Try Ensembl first
        if gene_id.upper().startswith('ENS'):
            # Validate Ensembl ID matches requested species
            ensembl_species_code = self._get_ensembl_species_from_id(gene_id)
            expected_species = self._get_species_prefix(species_normalized)
            
            if ensembl_species_code and expected_species:
                if ensembl_species_code != expected_species:
                    logger.warning(f"Ensembl ID {gene_id} appears to be for a different species "
                                 f"(ID prefix: {ensembl_species_code}, expected: {expected_species})")
                    logger.info(f"Attempting to look up gene symbol for {species_normalized} instead")
                    # Try to get the gene symbol from the Ensembl ID and search for it in the correct species
                    result = self._lookup_ensembl_and_find_ortholog(gene_id, species_normalized, expand_5p, expand_3p)
                else:
                    result = self._get_from_ensembl(gene_id, expand_5p, expand_3p)
            else:
                result = self._get_from_ensembl(gene_id, expand_5p, expand_3p)
        else:
            # Gene symbol search - uses species
            result = self._search_gene_by_symbol(gene_id, species_normalized)

        if result:
            # Add species to result
            result['requested_species'] = species_normalized
            # Cache the result
            try:
                os.makedirs(os.path.dirname(cache_file), exist_ok=True)
                with open(cache_file, 'w') as f:
                    json.dump(result, f, indent=2)
                logger.info(f"Cached sequence to: {cache_file}")
            except Exception as e:
                logger.warning(f"Could not cache sequence: {e}")
        else:
            logger.warning(f"No sequence found for {gene_id} in {species_normalized}")

        return result
    
    def _normalize_species_for_ensembl(self, species: str) -> str:
        """Normalize species name to Ensembl format."""
        species_map = {
            'human': 'homo_sapiens',
            'homo sapiens': 'homo_sapiens',
            'mouse': 'mus_musculus',
            'mus musculus': 'mus_musculus',
            'rat': 'rattus_norvegicus',
            'rattus norvegicus': 'rattus_norvegicus',
            'cat': 'felis_catus',
            'felis catus': 'felis_catus',
            'dog': 'canis_lupus_familiaris',
            'zebrafish': 'danio_rerio',
            'chicken': 'gallus_gallus',
            'pig': 'sus_scrofa',
            'cow': 'bos_taurus',
        }
        species_clean = species.lower().replace('_', ' ')
        return species_map.get(species_clean, species.lower().replace(' ', '_'))
    
    def _get_ensembl_species_from_id(self, ensembl_id: str) -> str:
        """Extract species code from Ensembl ID."""
        # Ensembl ID format: ENS[species code][feature type][11 digits]
        # Human: ENSG... (no species code), ENST...
        # Mouse: ENSMUSG..., ENSMUST...
        # Rat: ENSRNOG..., ENSRNOT...
        # Cat: ENSFCAG..., ENSFCAT...
        
        ensembl_id = ensembl_id.upper()
        if ensembl_id.startswith('ENSG') or ensembl_id.startswith('ENST'):
            return 'human'
        elif ensembl_id.startswith('ENSMUSG') or ensembl_id.startswith('ENSMUST'):
            return 'mouse'
        elif ensembl_id.startswith('ENSRNOG') or ensembl_id.startswith('ENSRNOT'):
            return 'rat'
        elif ensembl_id.startswith('ENSFCAG') or ensembl_id.startswith('ENSFCAT'):
            return 'cat'
        elif ensembl_id.startswith('ENSCAFG') or ensembl_id.startswith('ENSCAFT'):
            return 'dog'
        elif ensembl_id.startswith('ENSDARG') or ensembl_id.startswith('ENSDART'):
            return 'zebrafish'
        elif ensembl_id.startswith('ENSGALG') or ensembl_id.startswith('ENSGALT'):
            return 'chicken'
        elif ensembl_id.startswith('ENSSSCG') or ensembl_id.startswith('ENSSSCT'):
            return 'pig'
        return None
    
    def _get_species_prefix(self, species: str) -> str:
        """Get expected Ensembl ID prefix for species."""
        species_prefix_map = {
            'homo_sapiens': 'human',
            'mus_musculus': 'mouse',
            'rattus_norvegicus': 'rat',
            'felis_catus': 'cat',
            'canis_lupus_familiaris': 'dog',
            'danio_rerio': 'zebrafish',
            'gallus_gallus': 'chicken',
            'sus_scrofa': 'pig',
        }
        return species_prefix_map.get(species.lower(), None)
    
    def _lookup_ensembl_and_find_ortholog(self, ensembl_id: str, target_species: str,
                                          expand_5p: int, expand_3p: int) -> Optional[Dict]:
        """Look up gene from Ensembl ID and find ortholog in target species."""
        try:
            logger.info(f"Looking up gene symbol from {ensembl_id} to find in {target_species}")
            
            # First, get gene info from the Ensembl ID
            gene_url = f"{self.ENSEMBL_REST}/lookup/id/{ensembl_id}"
            response = requests.get(gene_url, headers={"Content-Type": "application/json"}, timeout=30)
            
            if response.status_code == 200:
                gene_info = response.json()
                gene_symbol = gene_info.get('display_name', '')
                
                if gene_symbol:
                    logger.info(f"Found gene symbol: {gene_symbol}, searching in {target_species}")
                    # Now search for this gene symbol in the target species
                    return self._search_gene_by_symbol(gene_symbol, target_species)
            
            logger.warning(f"Could not find ortholog for {ensembl_id} in {target_species}")
            return None
            
        except Exception as e:
            logger.error(f"Error looking up ortholog: {e}")
            return None
    
    def _get_from_ensembl(self, gene_id: str, expand_5p: int = 500, 
                         expand_3p: int = 500) -> Optional[Dict]:
        """Retrieve sequence from Ensembl REST API."""
        try:
            logger.info(f"Fetching gene details for {gene_id} from Ensembl...")
            
            # Get gene information
            gene_url = f"{self.ENSEMBL_REST}/lookup/id/{gene_id}?expand=1"
            headers = {"Content-Type": "application/json"}
            
            logger.debug(f"Gene lookup URL: {gene_url}")
            response = requests.get(gene_url, headers=headers, timeout=30)
            
            if response.status_code != 200:
                logger.error(f"Failed to get gene info: {response.status_code} - {response.text[:200]}")
                return None
                
            gene_info = response.json()
            logger.info(f"Got gene info: {gene_info.get('display_name', 'Unknown')}")
            
            # Calculate sequence region
            chrom = gene_info.get('seq_region_name', '')
            start = gene_info.get('start', 1)
            end = gene_info.get('end', start + 1000)
            species = gene_info.get('species', 'homo_sapiens')
            
            seq_start = max(1, start - expand_5p)
            seq_end = end + expand_3p
            
            logger.info(f"Fetching sequence region: {chrom}:{seq_start}-{seq_end}")
            
            # Get sequence - use correct endpoint format
            seq_url = f"{self.ENSEMBL_REST}/sequence/region/{species}/{chrom}:{seq_start}..{seq_end}:1"
            logger.debug(f"Sequence URL: {seq_url}")
            
            seq_response = requests.get(seq_url, headers={"Content-Type": "text/plain"}, timeout=60)
            
            if seq_response.status_code != 200:
                # Try alternative: get cDNA sequence from canonical transcript
                logger.warning(f"Region sequence failed ({seq_response.status_code}), trying cDNA...")
                transcripts = gene_info.get('Transcript', [])
                logger.info(f"Found {len(transcripts)} transcripts")
                
                canonical = next((t for t in transcripts if t.get('is_canonical')), 
                                transcripts[0] if transcripts else None)
                
                if canonical:
                    cdna_url = f"{self.ENSEMBL_REST}/sequence/id/{canonical['id']}?type=cdna"
                    logger.debug(f"cDNA URL: {cdna_url}")
                    cdna_response = requests.get(cdna_url, headers={"Content-Type": "text/plain"}, timeout=30)
                    
                    if cdna_response.status_code == 200:
                        sequence = cdna_response.text.strip().upper()
                        logger.info(f"Got cDNA sequence: {len(sequence)} bp")
                    else:
                        logger.error(f"Could not fetch cDNA: {cdna_response.status_code}")
                        return None
                else:
                    logger.error("No transcripts available")
                    return None
            else:
                sequence = seq_response.text.strip().upper()
                logger.info(f"Got genomic sequence: {len(sequence)} bp")
            
            return {
                'gene_id': gene_id,
                'gene_name': gene_info.get('display_name', gene_id),
                'gene_symbol': gene_info.get('display_name', gene_id),
                'ensembl_id': gene_id,
                'species': species,
                'chromosome': chrom,
                'start': start,
                'end': end,
                'strand': gene_info.get('strand', 1),
                'sequence': sequence,
                'sequence_start': seq_start,
                'sequence_end': seq_end,
                'assembly': gene_info.get('assembly_name', 'GRCh38'),
                'source': 'ensembl',
                'exons': self._extract_exons(gene_info)
            }
            
        except requests.exceptions.Timeout:
            logger.error(f"Timeout fetching from Ensembl")
            return None
        except Exception as e:
            logger.error(f"Error retrieving from Ensembl: {str(e)}", exc_info=True)
            return None
    
    def _search_gene_by_symbol(self, gene_symbol: str, species: str) -> Optional[Dict]:
        """Search for gene by symbol using Ensembl REST API and fetch sequence."""
        try:
            # Comprehensive species mapping - all Ensembl supported species
            species_map = {
                # Common names and scientific names
                'human': 'homo_sapiens',
                'homo_sapiens': 'homo_sapiens',
                'mouse': 'mus_musculus',
                'mus_musculus': 'mus_musculus',
                'rat': 'rattus_norvegicus',
                'rattus_norvegicus': 'rattus_norvegicus',
                'cat': 'felis_catus',
                'felis_catus': 'felis_catus',
                'dog': 'canis_lupus_familiaris',
                'canis_familiaris': 'canis_lupus_familiaris',
                'canis_lupus_familiaris': 'canis_lupus_familiaris',
                'zebrafish': 'danio_rerio',
                'danio_rerio': 'danio_rerio',
                'chicken': 'gallus_gallus',
                'gallus_gallus': 'gallus_gallus',
                'pig': 'sus_scrofa',
                'sus_scrofa': 'sus_scrofa',
                'cow': 'bos_taurus',
                'cattle': 'bos_taurus',
                'bos_taurus': 'bos_taurus',
                'sheep': 'ovis_aries',
                'ovis_aries': 'ovis_aries',
                'goat': 'capra_hircus',
                'capra_hircus': 'capra_hircus',
                'horse': 'equus_caballus',
                'equus_caballus': 'equus_caballus',
                'rabbit': 'oryctolagus_cuniculus',
                'oryctolagus_cuniculus': 'oryctolagus_cuniculus',
                'guinea_pig': 'cavia_porcellus',
                'cavia_porcellus': 'cavia_porcellus',
                'hamster': 'mesocricetus_auratus',
                'mesocricetus_auratus': 'mesocricetus_auratus',
                'frog': 'xenopus_tropicalis',
                'xenopus': 'xenopus_tropicalis',
                'xenopus_tropicalis': 'xenopus_tropicalis',
                'fruit_fly': 'drosophila_melanogaster',
                'drosophila': 'drosophila_melanogaster',
                'drosophila_melanogaster': 'drosophila_melanogaster',
                'worm': 'caenorhabditis_elegans',
                'c_elegans': 'caenorhabditis_elegans',
                'caenorhabditis_elegans': 'caenorhabditis_elegans',
                'yeast': 'saccharomyces_cerevisiae',
                'saccharomyces_cerevisiae': 'saccharomyces_cerevisiae',
                'macaque': 'macaca_mulatta',
                'rhesus': 'macaca_mulatta',
                'macaca_mulatta': 'macaca_mulatta',
                'cynomolgus': 'macaca_fascicularis',
                'macaca_fascicularis': 'macaca_fascicularis',
                'chimpanzee': 'pan_troglodytes',
                'pan_troglodytes': 'pan_troglodytes',
                'gorilla': 'gorilla_gorilla',
                'gorilla_gorilla': 'gorilla_gorilla',
                'orangutan': 'pongo_abelii',
                'pongo_abelii': 'pongo_abelii',
                'marmoset': 'callithrix_jacchus',
                'callithrix_jacchus': 'callithrix_jacchus',
                'bonobo': 'pan_paniscus',
                'pan_paniscus': 'pan_paniscus',
                'ferret': 'mustela_putorius_furo',
                'mustela_putorius_furo': 'mustela_putorius_furo',
                'medaka': 'oryzias_latipes',
                'oryzias_latipes': 'oryzias_latipes',
                'pufferfish': 'takifugu_rubripes',
                'fugu': 'takifugu_rubripes',
                'takifugu_rubripes': 'takifugu_rubripes',
                'stickleback': 'gasterosteus_aculeatus',
                'gasterosteus_aculeatus': 'gasterosteus_aculeatus',
                'salmon': 'salmo_salar',
                'salmo_salar': 'salmo_salar',
                'tilapia': 'oreochromis_niloticus',
                'oreochromis_niloticus': 'oreochromis_niloticus',
                'turkey': 'meleagris_gallopavo',
                'meleagris_gallopavo': 'meleagris_gallopavo',
                'duck': 'anas_platyrhynchos',
                'anas_platyrhynchos': 'anas_platyrhynchos',
                'platypus': 'ornithorhynchus_anatinus',
                'ornithorhynchus_anatinus': 'ornithorhynchus_anatinus',
                'opossum': 'monodelphis_domestica',
                'monodelphis_domestica': 'monodelphis_domestica',
                'elephant': 'loxodonta_africana',
                'loxodonta_africana': 'loxodonta_africana',
                'armadillo': 'dasypus_novemcinctus',
                'dasypus_novemcinctus': 'dasypus_novemcinctus',
                'alpaca': 'vicugna_pacos',
                'vicugna_pacos': 'vicugna_pacos',
                'dolphin': 'tursiops_truncatus',
                'tursiops_truncatus': 'tursiops_truncatus',
                'bat': 'myotis_lucifugus',
                'myotis_lucifugus': 'myotis_lucifugus',
                'shrew': 'sorex_araneus',
                'sorex_araneus': 'sorex_araneus',
                'hedgehog': 'erinaceus_europaeus',
                'erinaceus_europaeus': 'erinaceus_europaeus',
                'squirrel': 'ictidomys_tridecemlineatus',
                'ictidomys_tridecemlineatus': 'ictidomys_tridecemlineatus',
                'kangaroo_rat': 'dipodomys_ordii',
                'dipodomys_ordii': 'dipodomys_ordii',
                'pika': 'ochotona_princeps',
                'ochotona_princeps': 'ochotona_princeps',
                'tree_shrew': 'tupaia_belangeri',
                'tupaia_belangeri': 'tupaia_belangeri',
                'bushbaby': 'otolemur_garnettii',
                'otolemur_garnettii': 'otolemur_garnettii',
                'tarsier': 'tarsius_syrichta',
                'tarsius_syrichta': 'tarsius_syrichta',
                'gibbon': 'nomascus_leucogenys',
                'nomascus_leucogenys': 'nomascus_leucogenys',
                'sloth': 'choloepus_hoffmanni',
                'choloepus_hoffmanni': 'choloepus_hoffmanni',
                'lamprey': 'petromyzon_marinus',
                'petromyzon_marinus': 'petromyzon_marinus',
                'coelacanth': 'latimeria_chalumnae',
                'latimeria_chalumnae': 'latimeria_chalumnae',
                'ciona': 'ciona_intestinalis',
                'ciona_intestinalis': 'ciona_intestinalis',
                'sea_urchin': 'strongylocentrotus_purpuratus',
                'strongylocentrotus_purpuratus': 'strongylocentrotus_purpuratus',
            }

            ensembl_species = species_map.get(species.lower(), species.lower().replace(' ', '_'))

            # Try lookup/symbol endpoint first (more reliable)
            lookup_url = f"{self.ENSEMBL_REST}/lookup/symbol/{ensembl_species}/{gene_symbol}"
            headers = {"Content-Type": "application/json"}

            logger.info(f"Searching for {gene_symbol} in {ensembl_species}...")
            logger.debug(f"URL: {lookup_url}")

            response = requests.get(lookup_url, headers=headers, timeout=30)

            logger.info(f"Ensembl lookup response: {response.status_code}")

            if response.status_code == 200:
                gene_info = response.json()
                gene_id = gene_info.get('id', '')

                if gene_id:
                    logger.info(f"Found gene: {gene_id}")
                    # Now fetch the actual sequence using the gene ID
                    logger.info(f"Fetching sequence for {gene_id}...")
                    return self._get_from_ensembl(gene_id)

            # Fallback: Try xrefs/symbol endpoint
            logger.info("Trying xrefs/symbol endpoint...")
            search_url = f"{self.ENSEMBL_REST}/xrefs/symbol/{ensembl_species}/{gene_symbol}"
            response = requests.get(search_url, headers=headers, timeout=30)

            if response.status_code == 200:
                results = response.json()
                # Find the gene (not transcript)
                for result in results:
                    if result.get('type') == 'gene':
                        return self._get_from_ensembl(result['id'])

            logger.warning(f"Gene {gene_symbol} not found in Ensembl for {ensembl_species}")
            return None

        except requests.exceptions.Timeout:
            logger.error(f"Timeout searching for gene {gene_symbol}")
            return None
        except Exception as e:
            logger.error(f"Error searching gene: {str(e)}")
            return None
    
    def _extract_exons(self, gene_info: Dict) -> List[Dict]:
        """Extract exon information from gene data."""
        exons = []
        
        if 'Transcript' in gene_info:
            # Use canonical transcript or first transcript
            transcripts = gene_info['Transcript']
            if transcripts:
                transcript = transcripts[0]
                if 'Exon' in transcript:
                    for exon in transcript['Exon']:
                        exons.append({
                            'start': exon['start'],
                            'end': exon['end'],
                            'id': exon.get('id', '')
                        })
                        
        return exons
    
    def extract_exon_sequences(self, gene_data: Dict) -> List[Dict]:
        """
        Extract individual exon sequences from gene data.
        
        Args:
            gene_data: Gene data dictionary from get_gene_sequence
            
        Returns:
            List of exon dictionaries with sequences
        """
        exon_sequences = []
        full_sequence = gene_data['sequence']
        seq_start = gene_data['sequence_start']
        
        for exon in gene_data.get('exons', []):
            # Calculate relative positions
            relative_start = exon['start'] - seq_start
            relative_end = exon['end'] - seq_start + 1
            
            if 0 <= relative_start < len(full_sequence) and relative_end <= len(full_sequence):
                exon_seq = full_sequence[relative_start:relative_end]
                
                # Reverse complement if on negative strand
                if gene_data['strand'] == -1:
                    exon_seq = str(Seq(exon_seq).reverse_complement())
                    
                exon_sequences.append({
                    'exon_id': exon['id'],
                    'sequence': exon_seq,
                    'start': exon['start'],
                    'end': exon['end'],
                    'length': len(exon_seq)
                })
                
        return exon_sequences
    
    def find_target_regions(self, gene_data: Dict, target_type: str = 'exons') -> List[Tuple[int, int]]:
        """
        Find suitable target regions for gRNA design.
        
        Args:
            gene_data: Gene data dictionary
            target_type: Type of regions to target ('exons', 'all', 'custom')
            
        Returns:
            List of (start, end) tuples for target regions
        """
        regions = []
        
        if target_type == 'exons':
            # Target all exons
            for exon in gene_data.get('exons', []):
                regions.append((exon['start'], exon['end']))
                
        elif target_type == 'all':
            # Target entire gene region
            regions.append((gene_data['start'], gene_data['end']))
            
        return regions
    
    def get_sequence_context(self, sequence: str, position: int, 
                           upstream: int = 30, downstream: int = 30) -> str:
        """
        Extract sequence context around a position.
        
        Args:
            sequence: Full sequence
            position: Target position
            upstream: Bases to include upstream
            downstream: Bases to include downstream
            
        Returns:
            Context sequence
        """
        start = max(0, position - upstream)
        end = min(len(sequence), position + downstream)
        
        return sequence[start:end]
    
    def analyze_sequence_features(self, sequence: str) -> Dict:
        """
        Analyze sequence for relevant features.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Dictionary of sequence features
        """
        features = {
            'length': len(sequence),
            'gc_content': self._calculate_gc_content(sequence),
            'homopolymers': self._find_homopolymers(sequence),
            'repeats': self._find_repeats(sequence),
            'restriction_sites': self._find_restriction_sites(sequence)
        }
        
        return features
    
    def _calculate_gc_content(self, sequence: str) -> float:
        """Calculate GC content percentage."""
        gc_count = sequence.upper().count('G') + sequence.upper().count('C')
        return (gc_count / len(sequence)) * 100 if sequence else 0
    
    def _find_homopolymers(self, sequence: str, min_length: int = 4) -> List[Dict]:
        """Find homopolymer runs in sequence."""
        homopolymers = []
        sequence = sequence.upper()
        
        for base in 'ATCG':
            pattern = base * min_length
            pos = 0
            while True:
                pos = sequence.find(pattern, pos)
                if pos == -1:
                    break
                    
                # Extend to find full length
                length = min_length
                while pos + length < len(sequence) and sequence[pos + length] == base:
                    length += 1
                    
                homopolymers.append({
                    'base': base,
                    'position': pos,
                    'length': length
                })
                
                pos += length
                
        return homopolymers
    
    def _find_repeats(self, sequence: str, min_length: int = 6) -> List[Dict]:
        """Find repeat sequences."""
        repeats = []
        sequence = sequence.upper()
        
        # Simple repeat finder
        for length in range(min_length, min(20, len(sequence) // 2)):
            for i in range(len(sequence) - length):
                pattern = sequence[i:i+length]
                
                # Count occurrences
                count = 0
                pos = i
                positions = []
                
                while True:
                    pos = sequence.find(pattern, pos)
                    if pos == -1:
                        break
                    positions.append(pos)
                    pos += 1
                    count += 1
                    
                if count >= 2:
                    repeats.append({
                        'sequence': pattern,
                        'count': count,
                        'positions': positions
                    })
                    
        # Remove redundant repeats
        unique_repeats = []
        for repeat in repeats:
            is_redundant = False
            for other in unique_repeats:
                if repeat['sequence'] in other['sequence']:
                    is_redundant = True
                    break
            if not is_redundant:
                unique_repeats.append(repeat)
                
        return unique_repeats
    
    def _find_restriction_sites(self, sequence: str) -> List[Dict]:
        """Find common restriction enzyme sites."""
        # Common restriction sites
        restriction_sites = {
            'EcoRI': 'GAATTC',
            'BamHI': 'GGATCC',
            'HindIII': 'AAGCTT',
            'PstI': 'CTGCAG',
            'SalI': 'GTCGAC',
            'XbaI': 'TCTAGA',
            'NotI': 'GCGGCCGC',
            'XhoI': 'CTCGAG',
            'SacI': 'GAGCTC',
            'KpnI': 'GGTACC'
        }
        
        found_sites = []
        sequence = sequence.upper()
        
        for enzyme, site in restriction_sites.items():
            pos = 0
            while True:
                pos = sequence.find(site, pos)
                if pos == -1:
                    break
                    
                found_sites.append({
                    'enzyme': enzyme,
                    'site': site,
                    'position': pos
                })
                
                pos += 1
                
        return found_sites 