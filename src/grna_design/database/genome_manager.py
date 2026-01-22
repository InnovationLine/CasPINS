"""
Genome Manager Module
Handles genome data retrieval, caching, and indexing
"""

import os
import json
import gzip
import sqlite3
from typing import Dict, List, Optional, Tuple
import requests
from Bio import SeqIO
import hashlib
import pickle


class GenomeManager:
    """Manage genome sequences and provide fast access."""
    
    GENOME_SOURCES = {
        'ensembl': {
            'base_url': 'https://ftp.ensembl.org/pub',
            'species_map': {
                'human': {'name': 'homo_sapiens', 'assembly': 'GRCh38'},
                'mouse': {'name': 'mus_musculus', 'assembly': 'GRCm39'},
                'rat': {'name': 'rattus_norvegicus', 'assembly': 'mRatBN7.2'},
                'zebrafish': {'name': 'danio_rerio', 'assembly': 'GRCz11'}
            }
        },
        'ncbi': {
            'base_url': 'https://ftp.ncbi.nlm.nih.gov/genomes',
            'species_map': {
                'human': {'taxid': '9606', 'assembly': 'GCF_000001405.40'},
                'mouse': {'taxid': '10090', 'assembly': 'GCF_000001635.27'}
            }
        }
    }
    
    def __init__(self, cache_dir: str = "cache/genomes"):
        """
        Initialize genome manager.
        
        Args:
            cache_dir: Directory for caching genome data
        """
        self.cache_dir = cache_dir
        self.db_path = os.path.join(cache_dir, "genome_index.db")
        os.makedirs(cache_dir, exist_ok=True)
        
        # Initialize database
        self._init_database()
        
    def _init_database(self):
        """Initialize SQLite database for genome indexing."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        # Create tables
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS genomes (
                id INTEGER PRIMARY KEY,
                species TEXT,
                assembly TEXT,
                chromosome TEXT,
                length INTEGER,
                file_path TEXT,
                source TEXT,
                download_date TEXT,
                UNIQUE(species, assembly, chromosome)
            )
        """)
        
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS genome_features (
                id INTEGER PRIMARY KEY,
                genome_id INTEGER,
                feature_type TEXT,
                chromosome TEXT,
                start INTEGER,
                end INTEGER,
                strand TEXT,
                name TEXT,
                attributes TEXT,
                FOREIGN KEY (genome_id) REFERENCES genomes (id)
            )
        """)
        
        # Create indexes
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_features_position 
            ON genome_features (chromosome, start, end)
        """)
        
        conn.commit()
        conn.close()
    
    def get_genome(self, species: str, assembly: str = None) -> Optional[Dict]:
        """
        Get genome data for a species.
        
        Args:
            species: Species name
            assembly: Assembly version (optional)
            
        Returns:
            Dictionary with genome information
        """
        # Check if genome is cached
        genome_info = self._get_cached_genome(species, assembly)
        
        if not genome_info:
            # Download genome
            print(f"Downloading genome for {species}...")
            genome_info = self._download_genome(species, assembly)
            
        return genome_info
    
    def _get_cached_genome(self, species: str, assembly: str = None) -> Optional[Dict]:
        """Check if genome is cached."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        if assembly:
            cursor.execute("""
                SELECT * FROM genomes 
                WHERE species = ? AND assembly = ?
            """, (species, assembly))
        else:
            cursor.execute("""
                SELECT * FROM genomes 
                WHERE species = ? 
                ORDER BY download_date DESC 
                LIMIT 1
            """, (species,))
            
        result = cursor.fetchone()
        conn.close()
        
        if result:
            return {
                'species': result[1],
                'assembly': result[2],
                'chromosomes': self._get_genome_chromosomes(species, result[2])
            }
        
        return None
    
    def _get_genome_chromosomes(self, species: str, assembly: str) -> Dict:
        """Get all chromosomes for a genome."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        cursor.execute("""
            SELECT chromosome, file_path, length 
            FROM genomes 
            WHERE species = ? AND assembly = ?
        """, (species, assembly))
        
        chromosomes = {}
        for row in cursor.fetchall():
            chromosomes[row[0]] = {
                'file_path': row[1],
                'length': row[2]
            }
            
        conn.close()
        return chromosomes
    
    def _download_genome(self, species: str, assembly: str = None) -> Optional[Dict]:
        """Download genome from Ensembl or NCBI."""
        # For now, return None - in production, implement actual download
        # This would involve:
        # 1. Downloading FASTA files
        # 2. Indexing sequences
        # 3. Storing in database
        
        print(f"Note: Automatic genome download not implemented in this demo")
        print(f"Please manually download genome files for {species}")
        return None
    
    def load_custom_genome(self, fasta_file: str, species: str, 
                          assembly: str, source: str = "custom") -> bool:
        """
        Load a custom genome from FASTA file.
        
        Args:
            fasta_file: Path to FASTA file
            species: Species name
            assembly: Assembly name
            source: Source of genome
            
        Returns:
            Success status
        """
        try:
            import datetime
            
            # Parse FASTA file
            sequences = {}
            for record in SeqIO.parse(fasta_file, "fasta"):
                sequences[record.id] = str(record.seq)
                
            # Store sequences
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            for chrom, seq in sequences.items():
                # Save sequence to file
                seq_file = os.path.join(
                    self.cache_dir, 
                    f"{species}_{assembly}_{chrom}.pkl"
                )
                
                with open(seq_file, 'wb') as f:
                    pickle.dump(seq, f)
                    
                # Update database
                cursor.execute("""
                    INSERT OR REPLACE INTO genomes 
                    (species, assembly, chromosome, length, file_path, source, download_date)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                """, (species, assembly, chrom, len(seq), seq_file, 
                     source, datetime.datetime.now().isoformat()))
                     
            conn.commit()
            conn.close()
            
            print(f"Successfully loaded {len(sequences)} chromosomes for {species} {assembly}")
            return True
            
        except Exception as e:
            print(f"Error loading genome: {str(e)}")
            return False
    
    def get_sequence(self, species: str, chromosome: str, 
                    start: int = None, end: int = None, 
                    assembly: str = None) -> Optional[str]:
        """
        Get sequence for a genomic region.
        
        Args:
            species: Species name
            chromosome: Chromosome name
            start: Start position (1-based)
            end: End position (1-based)
            assembly: Assembly version
            
        Returns:
            DNA sequence or None
        """
        # Get genome info
        genome_info = self._get_cached_genome(species, assembly)
        if not genome_info:
            return None
            
        # Check if chromosome exists
        if chromosome not in genome_info['chromosomes']:
            # Try with 'chr' prefix
            if not chromosome.startswith('chr'):
                chromosome = 'chr' + chromosome
                if chromosome not in genome_info['chromosomes']:
                    return None
                    
        # Load sequence
        chrom_info = genome_info['chromosomes'][chromosome]
        
        try:
            with open(chrom_info['file_path'], 'rb') as f:
                sequence = pickle.load(f)
                
            # Extract region if specified
            if start is not None and end is not None:
                # Convert to 0-based
                start = max(0, start - 1)
                end = min(len(sequence), end)
                return sequence[start:end]
            else:
                return sequence
                
        except Exception as e:
            print(f"Error loading sequence: {str(e)}")
            return None
    
    def search_sequence(self, query: str, species: str, 
                       max_mismatches: int = 0, 
                       assembly: str = None) -> List[Dict]:
        """
        Search for sequence in genome.
        
        Args:
            query: Query sequence
            species: Species to search
            max_mismatches: Maximum allowed mismatches
            assembly: Assembly version
            
        Returns:
            List of matches
        """
        matches = []
        genome_info = self._get_cached_genome(species, assembly)
        
        if not genome_info:
            return matches
            
        # Search each chromosome
        for chrom, info in genome_info['chromosomes'].items():
            chrom_seq = self.get_sequence(species, chrom, assembly=assembly)
            if not chrom_seq:
                continue
                
            # Simple exact match search (extend for mismatches)
            pos = 0
            while True:
                pos = chrom_seq.find(query, pos)
                if pos == -1:
                    break
                    
                matches.append({
                    'chromosome': chrom,
                    'position': pos + 1,  # Convert to 1-based
                    'strand': '+',
                    'sequence': query
                })
                
                pos += 1
                
        return matches
    
    def get_gene_annotations(self, species: str, assembly: str = None) -> List[Dict]:
        """
        Get gene annotations for genome.
        
        Args:
            species: Species name
            assembly: Assembly version
            
        Returns:
            List of gene annotations
        """
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        # Get genome ID
        cursor.execute("""
            SELECT id FROM genomes 
            WHERE species = ? AND assembly = ? 
            LIMIT 1
        """, (species, assembly or ''))
        
        result = cursor.fetchone()
        if not result:
            conn.close()
            return []
            
        genome_id = result[0]
        
        # Get features
        cursor.execute("""
            SELECT * FROM genome_features 
            WHERE genome_id = ? AND feature_type = 'gene'
        """, (genome_id,))
        
        genes = []
        for row in cursor.fetchall():
            genes.append({
                'chromosome': row[3],
                'start': row[4],
                'end': row[5],
                'strand': row[6],
                'name': row[7],
                'attributes': json.loads(row[8]) if row[8] else {}
            })
            
        conn.close()
        return genes
    
    def add_genome_feature(self, species: str, assembly: str,
                          feature_type: str, chromosome: str,
                          start: int, end: int, strand: str,
                          name: str, attributes: Dict = None):
        """Add a genomic feature annotation."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        # Get genome ID
        cursor.execute("""
            SELECT id FROM genomes 
            WHERE species = ? AND assembly = ? AND chromosome = ?
        """, (species, assembly, chromosome))
        
        result = cursor.fetchone()
        if not result:
            conn.close()
            return False
            
        genome_id = result[0]
        
        # Insert feature
        cursor.execute("""
            INSERT INTO genome_features 
            (genome_id, feature_type, chromosome, start, end, strand, name, attributes)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """, (genome_id, feature_type, chromosome, start, end, strand, 
             name, json.dumps(attributes) if attributes else None))
             
        conn.commit()
        conn.close()
        return True
    
    def get_genome_stats(self) -> Dict:
        """Get statistics about cached genomes."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        # Count genomes
        cursor.execute("SELECT COUNT(DISTINCT species || '_' || assembly) FROM genomes")
        num_genomes = cursor.fetchone()[0]
        
        # Get species list
        cursor.execute("SELECT DISTINCT species, assembly FROM genomes")
        species_list = cursor.fetchall()
        
        # Calculate total size
        total_size = 0
        for species, assembly in species_list:
            cursor.execute("""
                SELECT SUM(length) FROM genomes 
                WHERE species = ? AND assembly = ?
            """, (species, assembly))
            size = cursor.fetchone()[0]
            if size:
                total_size += size
                
        conn.close()
        
        return {
            'num_genomes': num_genomes,
            'species': [{'species': s[0], 'assembly': s[1]} for s in species_list],
            'total_size': total_size,
            'cache_dir': self.cache_dir
        } 