"""
Configuration module for CRISPR Analysis Suite

Provides:
- Species configuration (90+ species)
- CRISPR system configurations
- Editing mode configurations
"""

import json
import os
from typing import Dict, List, Optional

# Path to configuration files
CONFIG_DIR = os.path.dirname(os.path.abspath(__file__))
SPECIES_CONFIG_FILE = os.path.join(CONFIG_DIR, 'species.json')

# Cache for loaded configurations
_species_config_cache = None


def load_species_config() -> Dict:
    """
    Load the comprehensive species configuration.
    
    Returns:
        Dictionary containing all species, assemblies, and CRISPR systems
    """
    global _species_config_cache
    
    if _species_config_cache is not None:
        return _species_config_cache
    
    if os.path.exists(SPECIES_CONFIG_FILE):
        with open(SPECIES_CONFIG_FILE, 'r') as f:
            _species_config_cache = json.load(f)
    else:
        # Fallback minimal config
        _species_config_cache = {
            "species": {
                "mammals": [
                    {
                        "common_name": "Human",
                        "scientific_name": "Homo sapiens",
                        "ensembl_name": "homo_sapiens",
                        "assemblies": [{"name": "GRCh38", "alias": "hg38", "default": True}]
                    }
                ]
            },
            "crispr_systems": {
                "cas9_variants": [
                    {"name": "SpCas9", "pam": "NGG", "default": True}
                ]
            }
        }
    
    return _species_config_cache


def get_all_species() -> List[Dict]:
    """
    Get flat list of all supported species.
    
    Returns:
        List of species dictionaries with display names
    """
    config = load_species_config()
    all_species = []
    
    for category, species_list in config.get('species', {}).items():
        for species in species_list:
            display_name = f"{species['common_name']} ({species['scientific_name']})"
            default_assembly = next(
                (a['name'] for a in species.get('assemblies', []) if a.get('default')),
                species.get('assemblies', [{}])[0].get('name', 'Unknown')
            )
            
            all_species.append({
                'display_name': display_name,
                'common_name': species['common_name'],
                'scientific_name': species['scientific_name'],
                'ensembl_name': species.get('ensembl_name', ''),
                'category': category,
                'default_assembly': default_assembly,
                'assemblies': species.get('assemblies', []),
                'ncbi_taxonomy_id': species.get('ncbi_taxonomy_id')
            })
    
    return all_species


def get_species_for_dropdown() -> List[str]:
    """
    Get list of species names formatted for dropdown selection.
    
    Returns:
        List of formatted species strings
    """
    all_species = get_all_species()
    dropdown_options = []
    
    # Group by category for organized dropdown
    categories = {}
    for sp in all_species:
        cat = sp['category'].replace('_', ' ').title()
        if cat not in categories:
            categories[cat] = []
        categories[cat].append(sp['display_name'])
    
    # Add category headers and species
    for category in ['Mammals', 'Birds', 'Fish', 'Amphibians', 'Reptiles', 
                     'Invertebrates', 'Plants', 'Fungi', 'Bacteria']:
        if category in categories:
            for species in sorted(categories[category]):
                dropdown_options.append(species)
    
    return dropdown_options


def get_crispr_systems() -> List[Dict]:
    """
    Get list of all supported CRISPR systems.
    
    Returns:
        List of CRISPR system dictionaries
    """
    config = load_species_config()
    systems = []
    
    for category, system_list in config.get('crispr_systems', {}).items():
        for system in system_list:
            systems.append({
                'name': system['name'],
                'pam': system['pam'],
                'default': system.get('default', False),
                'category': category.replace('_', ' ').title()
            })
    
    return systems


def get_editing_modes() -> List[Dict]:
    """
    Get list of available editing modes.
    
    Returns:
        List of editing mode dictionaries
    """
    config = load_species_config()
    return config.get('editing_modes', [
        {"name": "knockout", "description": "Gene disruption via NHEJ", "default": True}
    ])


def get_species_by_name(name: str) -> Optional[Dict]:
    """
    Find species by common or scientific name.
    
    Args:
        name: Species name to search for
        
    Returns:
        Species dictionary or None if not found
    """
    all_species = get_all_species()
    name_lower = name.lower()
    
    for sp in all_species:
        if (sp['common_name'].lower() == name_lower or 
            sp['scientific_name'].lower() == name_lower or
            sp['ensembl_name'].lower() == name_lower):
            return sp
    
    # Partial match
    for sp in all_species:
        if (name_lower in sp['common_name'].lower() or 
            name_lower in sp['scientific_name'].lower()):
            return sp
    
    return None


def normalize_species_name(name: str) -> str:
    """
    Convert various species name formats to ensembl name.
    
    Args:
        name: Input species name (e.g., "Human", "Homo sapiens", "human")
        
    Returns:
        Ensembl-compatible species name (e.g., "homo_sapiens")
    """
    species = get_species_by_name(name)
    if species:
        return species['ensembl_name']
    
    # Fallback: basic normalization
    return name.lower().replace(' ', '_')


# Import settings functions
from .settings import (
    get_data_directory,
    set_data_directory,
    is_data_directory_configured,
    list_available_genes,
    get_gene_folder,
    load_config,
    save_config,
    get_cache_directory,
    get_output_directory
)


# Convenience exports
__all__ = [
    # Species config
    'load_species_config',
    'get_all_species',
    'get_species_for_dropdown',
    'get_crispr_systems',
    'get_editing_modes',
    'get_species_by_name',
    'normalize_species_name',
    # Settings/paths
    'get_data_directory',
    'set_data_directory',
    'is_data_directory_configured',
    'list_available_genes',
    'get_gene_folder',
    'load_config',
    'save_config',
    'get_cache_directory',
    'get_output_directory'
]
