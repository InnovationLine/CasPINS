"""
HDR (Homology-Directed Repair) Template Designer Module

Designs donor templates for knock-in experiments using CRISPR-Cas9.
Includes:
- Homology arm design (5' and 3')
- Silent mutation insertion to prevent re-cutting
- ssODN and dsDNA template generation
- Integration with gRNA design for optimal cut site selection
"""

import re
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass
from enum import Enum


class HDRTemplateType(Enum):
    """Types of HDR donor templates."""
    SSODN = "ssODN"  # Single-stranded oligodeoxynucleotide (small edits)
    DSDNA = "dsDNA"  # Double-stranded DNA (larger insertions)
    PLASMID = "plasmid"  # Plasmid donor (large insertions, selection markers)


class EditType(Enum):
    """Types of edits for HDR."""
    POINT_MUTATION = "point_mutation"
    SMALL_INSERTION = "small_insertion"  # < 50bp
    SMALL_DELETION = "small_deletion"    # < 50bp
    TAG_INSERTION = "tag_insertion"      # GFP, FLAG, HA, etc.
    GENE_REPLACEMENT = "gene_replacement"
    KNOCK_IN = "knock_in"


@dataclass
class HDRTemplate:
    """Data class for HDR template information."""
    template_type: HDRTemplateType
    edit_type: EditType
    left_homology_arm: str
    right_homology_arm: str
    insert_sequence: str
    full_template: str
    target_sequence: str
    cut_site: int
    silent_mutations: List[Dict]
    warnings: List[str]
    
    # Metadata
    left_arm_length: int = 0
    right_arm_length: int = 0
    total_length: int = 0
    gc_content: float = 0.0


class HDRDesigner:
    """
    Designer for HDR (Homology-Directed Repair) templates.
    
    Supports:
    - Point mutations
    - Small insertions/deletions
    - Tag insertions (GFP, FLAG, HA, V5, Myc, His)
    - Gene knock-ins
    """
    
    # Common protein tags
    PROTEIN_TAGS = {
        'FLAG': 'GACTACAAAGACGATGACGACAAG',  # DYKDDDDK
        'HA': 'TACCCATACGATGTTCCAGATTACGCT',  # YPYDVPDYA
        'Myc': 'GAACAAAAACTCATCTCAGAAGAGGATCTG',  # EQKLISEEDL
        'V5': 'GGTAAGCCTATCCCTAACCCTCTCCTCGGTCTCGATTCTACG',  # GKPIPNPLLGLDST
        'His6': 'CATCATCACCATCACCAC',  # HHHHHH
        'His10': 'CATCACCATCACCATCACCATCACCATCAC',  # HHHHHHHHHH
        'GFP': None,  # Too large, requires plasmid
        'mCherry': None,  # Too large, requires plasmid
        'T2A': 'GAGGGCAGGGGAAGTCTTCTAACATGCGGGGACGTGGAGGAAAATCCCGGCCCC',  # Self-cleaving peptide
        'P2A': 'GCCACAAACTTCTCTCTGCTAAAGCAAGCAGGTGATGTTGAAGAAAACCCCGGGCCC',
        'GSG_linker': 'GGTTCTGGT',  # Flexible linker
        'GGGGS_linker': 'GGTGGAGGCGGTAGCGGTGGAGGCGGTAGC',  # Longer flexible linker
    }
    
    # Codon table for silent mutations
    CODON_TABLE = {
        'F': ['TTT', 'TTC'],
        'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
        'I': ['ATT', 'ATC', 'ATA'],
        'M': ['ATG'],
        'V': ['GTT', 'GTC', 'GTA', 'GTG'],
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
        'P': ['CCT', 'CCC', 'CCA', 'CCG'],
        'T': ['ACT', 'ACC', 'ACA', 'ACG'],
        'A': ['GCT', 'GCC', 'GCA', 'GCG'],
        'Y': ['TAT', 'TAC'],
        '*': ['TAA', 'TAG', 'TGA'],
        'H': ['CAT', 'CAC'],
        'Q': ['CAA', 'CAG'],
        'N': ['AAT', 'AAC'],
        'K': ['AAA', 'AAG'],
        'D': ['GAT', 'GAC'],
        'E': ['GAA', 'GAG'],
        'C': ['TGT', 'TGC'],
        'W': ['TGG'],
        'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    }
    
    # Reverse codon table
    REVERSE_CODON = {}
    for aa, codons in CODON_TABLE.items():
        for codon in codons:
            REVERSE_CODON[codon] = aa
    
    def __init__(self, species: str = 'human'):
        """
        Initialize HDR Designer.
        
        Args:
            species: Target species for codon optimization
        """
        self.species = species
        
    def design_hdr_template(self,
                           target_sequence: str,
                           cut_site: int,
                           edit_type: EditType,
                           insert_sequence: str = "",
                           left_arm_length: int = 800,
                           right_arm_length: int = 800,
                           add_silent_mutations: bool = True,
                           grna_sequence: str = None,
                           pam_sequence: str = "NGG") -> HDRTemplate:
        """
        Design an HDR donor template.
        
        Args:
            target_sequence: Full target genomic sequence
            cut_site: Position of Cas9 cut site in target sequence
            edit_type: Type of edit to make
            insert_sequence: Sequence to insert (for insertions)
            left_arm_length: Length of left homology arm (bp)
            right_arm_length: Length of right homology arm (bp)
            add_silent_mutations: Whether to add silent mutations to prevent re-cutting
            grna_sequence: gRNA sequence (for silent mutation design)
            pam_sequence: PAM sequence
            
        Returns:
            HDRTemplate object with complete donor design
        """
        warnings = []
        
        # Validate inputs
        if cut_site < 0 or cut_site > len(target_sequence):
            raise ValueError(f"Cut site {cut_site} is outside target sequence")
        
        # Determine template type based on edit
        if edit_type == EditType.POINT_MUTATION:
            template_type = HDRTemplateType.SSODN
            if left_arm_length > 60:
                left_arm_length = 60
                right_arm_length = 60
                warnings.append("Reduced arm length to 60bp for ssODN template")
        elif edit_type in [EditType.SMALL_INSERTION, EditType.SMALL_DELETION]:
            template_type = HDRTemplateType.SSODN
            if left_arm_length > 80:
                left_arm_length = 80
                right_arm_length = 80
        elif edit_type == EditType.TAG_INSERTION:
            # Check if tag is small enough for ssODN
            if len(insert_sequence) < 100:
                template_type = HDRTemplateType.SSODN
                left_arm_length = min(80, left_arm_length)
                right_arm_length = min(80, right_arm_length)
            else:
                template_type = HDRTemplateType.DSDNA
        else:
            template_type = HDRTemplateType.DSDNA
            if len(insert_sequence) > 2000:
                template_type = HDRTemplateType.PLASMID
                warnings.append("Large insert - consider using plasmid donor with selection marker")
        
        # Extract homology arms
        left_start = max(0, cut_site - left_arm_length)
        left_arm = target_sequence[left_start:cut_site]
        
        right_end = min(len(target_sequence), cut_site + right_arm_length)
        right_arm = target_sequence[cut_site:right_end]
        
        # Adjust arm lengths if near sequence boundaries
        actual_left_length = len(left_arm)
        actual_right_length = len(right_arm)
        
        if actual_left_length < left_arm_length:
            warnings.append(f"Left arm truncated to {actual_left_length}bp (near sequence start)")
        if actual_right_length < right_arm_length:
            warnings.append(f"Right arm truncated to {actual_right_length}bp (near sequence end)")
        
        # Design silent mutations if requested
        silent_mutations = []
        if add_silent_mutations and grna_sequence:
            left_arm, right_arm, silent_mutations = self._add_silent_mutations(
                left_arm, right_arm, grna_sequence, pam_sequence, cut_site
            )
            if silent_mutations:
                warnings.append(f"Added {len(silent_mutations)} silent mutation(s) to prevent re-cutting")
        
        # Construct full template
        full_template = left_arm + insert_sequence + right_arm
        
        # Calculate GC content
        gc_count = full_template.count('G') + full_template.count('C')
        gc_content = (gc_count / len(full_template) * 100) if full_template else 0
        
        # Check for issues
        if gc_content > 65:
            warnings.append(f"High GC content ({gc_content:.1f}%) may affect synthesis")
        elif gc_content < 35:
            warnings.append(f"Low GC content ({gc_content:.1f}%) may affect stability")
        
        # Check for problematic sequences
        if 'GGGG' in full_template or 'CCCC' in full_template:
            warnings.append("Homopolymer runs detected - may affect synthesis")
        
        return HDRTemplate(
            template_type=template_type,
            edit_type=edit_type,
            left_homology_arm=left_arm,
            right_homology_arm=right_arm,
            insert_sequence=insert_sequence,
            full_template=full_template,
            target_sequence=target_sequence[left_start:right_end],
            cut_site=cut_site,
            silent_mutations=silent_mutations,
            warnings=warnings,
            left_arm_length=actual_left_length,
            right_arm_length=actual_right_length,
            total_length=len(full_template),
            gc_content=gc_content
        )
    
    def design_tag_insertion(self,
                            target_sequence: str,
                            cut_site: int,
                            tag_name: str,
                            position: str = "C-terminal",
                            add_linker: bool = True,
                            grna_sequence: str = None) -> HDRTemplate:
        """
        Design HDR template for protein tag insertion.
        
        Args:
            target_sequence: Target genomic sequence
            cut_site: Position to insert tag
            tag_name: Name of tag (FLAG, HA, Myc, V5, His6, etc.)
            position: "N-terminal", "C-terminal", or "internal"
            add_linker: Whether to add flexible linker
            grna_sequence: gRNA sequence for silent mutations
            
        Returns:
            HDRTemplate for tag insertion
        """
        # Get tag sequence
        tag_name_upper = tag_name.upper()
        if tag_name_upper not in self.PROTEIN_TAGS:
            raise ValueError(f"Unknown tag: {tag_name}. Available: {list(self.PROTEIN_TAGS.keys())}")
        
        tag_seq = self.PROTEIN_TAGS[tag_name_upper]
        if tag_seq is None:
            raise ValueError(f"{tag_name} is too large for ssODN. Use plasmid donor instead.")
        
        # Add linker if requested
        insert_seq = ""
        if add_linker:
            linker = self.PROTEIN_TAGS['GSG_linker']
            if position == "N-terminal":
                insert_seq = tag_seq + linker
            elif position == "C-terminal":
                insert_seq = linker + tag_seq
            else:
                insert_seq = linker + tag_seq + linker
        else:
            insert_seq = tag_seq
        
        return self.design_hdr_template(
            target_sequence=target_sequence,
            cut_site=cut_site,
            edit_type=EditType.TAG_INSERTION,
            insert_sequence=insert_seq,
            left_arm_length=80,
            right_arm_length=80,
            grna_sequence=grna_sequence
        )
    
    def design_point_mutation(self,
                             target_sequence: str,
                             mutation_position: int,
                             original_base: str,
                             new_base: str,
                             grna_sequence: str = None) -> HDRTemplate:
        """
        Design ssODN for point mutation.
        
        Args:
            target_sequence: Target sequence
            mutation_position: Position of mutation
            original_base: Original nucleotide
            new_base: New nucleotide
            grna_sequence: gRNA sequence
            
        Returns:
            HDRTemplate for point mutation
        """
        # Verify original base
        if target_sequence[mutation_position] != original_base:
            raise ValueError(f"Expected {original_base} at position {mutation_position}, "
                           f"found {target_sequence[mutation_position]}")
        
        # Create modified sequence
        modified_seq = (target_sequence[:mutation_position] + 
                       new_base + 
                       target_sequence[mutation_position + 1:])
        
        # Design template centered on mutation
        arm_length = 60
        cut_site = mutation_position  # Approximate
        
        left_start = max(0, mutation_position - arm_length)
        right_end = min(len(target_sequence), mutation_position + arm_length + 1)
        
        left_arm = modified_seq[left_start:mutation_position]
        right_arm = modified_seq[mutation_position + 1:right_end]
        
        return HDRTemplate(
            template_type=HDRTemplateType.SSODN,
            edit_type=EditType.POINT_MUTATION,
            left_homology_arm=left_arm,
            right_homology_arm=right_arm,
            insert_sequence=new_base,
            full_template=left_arm + new_base + right_arm,
            target_sequence=target_sequence[left_start:right_end],
            cut_site=mutation_position,
            silent_mutations=[],
            warnings=[],
            left_arm_length=len(left_arm),
            right_arm_length=len(right_arm),
            total_length=len(left_arm) + 1 + len(right_arm),
            gc_content=0.0  # Calculate if needed
        )
    
    def _add_silent_mutations(self,
                             left_arm: str,
                             right_arm: str,
                             grna_sequence: str,
                             pam_sequence: str,
                             cut_site: int) -> Tuple[str, str, List[Dict]]:
        """
        Add silent mutations to PAM or seed region to prevent re-cutting.
        
        Strategy:
        1. First try to mutate PAM (NGG -> NAG, NCG, etc.)
        2. If PAM can't be mutated, add silent mutations to seed region
        """
        mutations = []
        
        # Find gRNA binding site in homology arms
        combined = left_arm + right_arm
        grna_pos = combined.find(grna_sequence)
        
        if grna_pos == -1:
            # Try reverse complement
            rev_comp = self._reverse_complement(grna_sequence)
            grna_pos = combined.find(rev_comp)
        
        if grna_pos == -1:
            return left_arm, right_arm, []  # gRNA not found in arms
        
        # PAM is typically 3bp after gRNA (for SpCas9 NGG)
        pam_pos = grna_pos + len(grna_sequence)
        
        # Try to mutate PAM
        if pam_pos + 3 <= len(combined):
            pam = combined[pam_pos:pam_pos + 3]
            
            # For NGG PAM, mutate second G to A, C, or T
            if pam[1:] == 'GG':
                # G -> A mutation in PAM creates NAG (much less active)
                new_pam = pam[0] + 'AG'
                combined = combined[:pam_pos] + new_pam + combined[pam_pos + 3:]
                mutations.append({
                    'position': pam_pos + 1,
                    'original': 'G',
                    'mutated': 'A',
                    'type': 'PAM_disruption',
                    'note': 'NGG -> NAG (>90% activity reduction)'
                })
        
        # Split back into arms
        left_arm = combined[:len(left_arm)]
        right_arm = combined[len(left_arm):]
        
        return left_arm, right_arm, mutations
    
    def _reverse_complement(self, sequence: str) -> str:
        """Get reverse complement of DNA sequence."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(base, 'N') for base in reversed(sequence))
    
    def get_ordering_info(self, template: HDRTemplate) -> Dict:
        """
        Get information for ordering the HDR template.
        
        Args:
            template: HDRTemplate object
            
        Returns:
            Dictionary with ordering information
        """
        info = {
            'template_type': template.template_type.value,
            'sequence': template.full_template,
            'length': template.total_length,
            'gc_content': f"{template.gc_content:.1f}%",
            'modifications': [],
            'notes': []
        }
        
        if template.template_type == HDRTemplateType.SSODN:
            info['notes'].append("Order as ultramer/megamer oligo")
            info['notes'].append("Consider phosphorothioate modifications at ends for stability")
            info['modifications'].append("5' and 3' phosphorothioate bonds (optional)")
            
            if template.total_length > 200:
                info['notes'].append("Length >200nt - verify synthesis capability with vendor")
        
        elif template.template_type == HDRTemplateType.DSDNA:
            info['notes'].append("Order as gene synthesis fragment")
            info['notes'].append("Can be PCR amplified from plasmid")
        
        elif template.template_type == HDRTemplateType.PLASMID:
            info['notes'].append("Clone into appropriate vector")
            info['notes'].append("Consider adding selection marker for enrichment")
        
        return info
    
    def format_for_display(self, template: HDRTemplate) -> str:
        """
        Format HDR template for display.
        
        Args:
            template: HDRTemplate object
            
        Returns:
            Formatted string for display
        """
        lines = []
        lines.append("=" * 80)
        lines.append(f"HDR DONOR TEMPLATE - {template.edit_type.value.upper()}")
        lines.append("=" * 80)
        lines.append("")
        
        lines.append(f"Template Type: {template.template_type.value}")
        lines.append(f"Total Length: {template.total_length} bp")
        lines.append(f"GC Content: {template.gc_content:.1f}%")
        lines.append("")
        
        lines.append("HOMOLOGY ARMS:")
        lines.append(f"  Left arm: {template.left_arm_length} bp")
        lines.append(f"  Right arm: {template.right_arm_length} bp")
        lines.append("")
        
        if template.insert_sequence:
            lines.append(f"INSERT SEQUENCE ({len(template.insert_sequence)} bp):")
            # Display insert in chunks of 60
            for i in range(0, len(template.insert_sequence), 60):
                lines.append(f"  {template.insert_sequence[i:i+60]}")
            lines.append("")
        
        lines.append("FULL TEMPLATE SEQUENCE:")
        lines.append("-" * 80)
        
        # Format sequence with annotations
        seq = template.full_template
        left_len = template.left_arm_length
        insert_len = len(template.insert_sequence)
        
        # Add position markers
        for i in range(0, len(seq), 60):
            chunk = seq[i:i+60]
            pos_str = str(i + 1).rjust(6)
            lines.append(f"{pos_str}  {chunk}")
        
        lines.append("-" * 80)
        lines.append("")
        
        # Annotations
        lines.append("SEQUENCE MAP:")
        lines.append(f"  [1-{left_len}] = Left homology arm")
        if insert_len > 0:
            lines.append(f"  [{left_len + 1}-{left_len + insert_len}] = Insert")
        lines.append(f"  [{left_len + insert_len + 1}-{template.total_length}] = Right homology arm")
        lines.append("")
        
        if template.silent_mutations:
            lines.append("SILENT MUTATIONS (to prevent re-cutting):")
            for mut in template.silent_mutations:
                lines.append(f"  Position {mut['position']}: {mut['original']} -> {mut['mutated']} ({mut['type']})")
            lines.append("")
        
        if template.warnings:
            lines.append("WARNINGS:")
            for warning in template.warnings:
                lines.append(f"  ⚠ {warning}")
            lines.append("")
        
        return "\n".join(lines)


# Convenience function for quick template design
def design_knockin_template(target_sequence: str,
                           cut_site: int,
                           insert_sequence: str,
                           grna_sequence: str = None,
                           species: str = 'human') -> HDRTemplate:
    """
    Quick function to design a knock-in HDR template.
    
    Args:
        target_sequence: Target genomic sequence
        cut_site: Cas9 cut site position
        insert_sequence: Sequence to insert
        grna_sequence: gRNA sequence (for silent mutations)
        species: Target species
        
    Returns:
        HDRTemplate object
    """
    designer = HDRDesigner(species)
    return designer.design_hdr_template(
        target_sequence=target_sequence,
        cut_site=cut_site,
        edit_type=EditType.KNOCK_IN,
        insert_sequence=insert_sequence,
        grna_sequence=grna_sequence
    )


__all__ = [
    'HDRDesigner',
    'HDRTemplate',
    'HDRTemplateType',
    'EditType',
    'design_knockin_template'
]
