"""
File Management Module
Handles archiving output files and directory operations
"""

import os
import shutil
from datetime import datetime


def archive_output_files(gene_folder):
    """
    Move all files from output folder to archive folder.
    Creates archive folder with timestamp if needed.
    
    Args:
        gene_folder: Path to the gene folder
        
    Returns:
        str: Archive directory path or None if nothing was archived
    """
    output_dir = os.path.join(gene_folder, "output")
    if not os.path.exists(output_dir):
        return None  # Nothing to archive
    
    # Get list of files in output directory
    files_to_archive = [f for f in os.listdir(output_dir) if os.path.isfile(os.path.join(output_dir, f))]
    
    if not files_to_archive:
        return None  # No files to archive
    
    # Create archive directory with timestamp
    archive_base = os.path.join(gene_folder, "archive")
    if not os.path.exists(archive_base):
        os.makedirs(archive_base)
    
    # Create timestamped subfolder in archive
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    archive_dir = os.path.join(archive_base, f"run_{timestamp}")
    os.makedirs(archive_dir)
    
    # Move all files from output to archive
    print(f"  Archiving {len(files_to_archive)} files from previous run...")
    for filename in files_to_archive:
        src = os.path.join(output_dir, filename)
        dst = os.path.join(archive_dir, filename)
        shutil.move(src, dst)
        print(f"    Archived: {filename}")
    
    print(f"  Files archived to: archive/run_{timestamp}/")
    return archive_dir


def ensure_output_dir(gene_folder):
    """
    Ensure output directory exists.
    
    Args:
        gene_folder: Path to the gene folder
        
    Returns:
        str: Path to output directory
    """
    output_dir = os.path.join(gene_folder, "output")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    return output_dir


def save_results(output_dir, filename, content, mode='w'):
    """
    Save results to a file.
    
    Args:
        output_dir: Output directory path
        filename: Name of the file
        content: Content to write
        mode: File open mode (default 'w')
        
    Returns:
        str: Full path to saved file
    """
    file_path = os.path.join(output_dir, filename)
    with open(file_path, mode, encoding='utf-8') as f:
        f.write(content)
    return file_path


def check_required_files(gene_folder):
    """
    Check if all required files exist in gene folder.
    
    Args:
        gene_folder: Path to the gene folder
        
    Returns:
        tuple: (bool, list) - (all_exist, missing_files)
    """
    required_files = {
        "control.ab1": "Control AB1 file",
        "edited.ab1": "Edited AB1 file", 
        "grna.txt": "gRNA sequences",
        "mrna.txt": "mRNA reference sequence"
    }
    
    missing_files = []
    for filename, description in required_files.items():
        file_path = os.path.join(gene_folder, filename)
        if not os.path.exists(file_path):
            missing_files.append(f"{filename} ({description})")
    
    return len(missing_files) == 0, missing_files 