"""
Application Settings and Configuration
Manages data paths and other configurable settings
"""

import os
import json
from pathlib import Path
from typing import Optional


# Default configuration
DEFAULT_CONFIG = {
    "data_directory": "",  # Empty means user must provide path
    "cache_directory": "cache",
    "output_directory": "output",
    "default_species": "human",
    "default_assembly": "GRCh38"
}

# Config file path
CONFIG_FILE = "crispr_config.json"


def get_config_path() -> Path:
    """Get path to config file."""
    # Check for config in current directory first
    local_config = Path(CONFIG_FILE)
    if local_config.exists():
        return local_config
    
    # Check in user's home directory
    home_config = Path.home() / ".crispr_analysis" / CONFIG_FILE
    return home_config


def load_config() -> dict:
    """Load configuration from file or return defaults."""
    config_path = get_config_path()
    
    if config_path.exists():
        try:
            with open(config_path, 'r') as f:
                config = json.load(f)
                # Merge with defaults to ensure all keys exist
                return {**DEFAULT_CONFIG, **config}
        except Exception as e:
            print(f"Warning: Could not load config: {e}")
    
    return DEFAULT_CONFIG.copy()


def save_config(config: dict):
    """Save configuration to file."""
    config_path = get_config_path()
    
    # Ensure directory exists
    config_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(config_path, 'w') as f:
        json.dump(config, f, indent=2)


def get_data_directory() -> Optional[str]:
    """Get the configured data directory."""
    config = load_config()
    data_dir = config.get("data_directory", "")
    
    if data_dir and os.path.isdir(data_dir):
        return data_dir
    
    return None


def set_data_directory(path: str) -> bool:
    """Set the data directory path."""
    if os.path.isdir(path):
        config = load_config()
        config["data_directory"] = path
        save_config(config)
        return True
    return False


def get_cache_directory() -> str:
    """Get cache directory path."""
    config = load_config()
    return config.get("cache_directory", "cache")


def get_output_directory() -> str:
    """Get output directory path."""
    config = load_config()
    return config.get("output_directory", "output")


def get_gene_folder(gene_name: str) -> Optional[str]:
    """Get path to a gene's data folder."""
    data_dir = get_data_directory()
    if not data_dir:
        return None
    
    gene_folder = os.path.join(data_dir, gene_name.lower())
    if os.path.isdir(gene_folder):
        return gene_folder
    
    return None


def list_available_genes() -> list:
    """List all genes available in the data directory."""
    data_dir = get_data_directory()
    if not data_dir:
        return []
    
    genes = []
    for item in os.listdir(data_dir):
        item_path = os.path.join(data_dir, item)
        if os.path.isdir(item_path):
            # Check if it has expected files
            has_grna = os.path.exists(os.path.join(item_path, "grna.txt"))
            has_mrna = os.path.exists(os.path.join(item_path, "mrna.txt"))
            has_input = os.path.isdir(os.path.join(item_path, "input"))
            
            genes.append({
                'name': item.upper(),
                'folder': item_path,
                'has_grna': has_grna,
                'has_mrna': has_mrna,
                'has_input': has_input
            })
    
    return genes


def is_data_directory_configured() -> bool:
    """Check if data directory is properly configured."""
    return get_data_directory() is not None
