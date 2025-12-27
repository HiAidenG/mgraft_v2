"""
CRISPR Clustering Module

Clusters CRISPR arrays based on:
- Repeat sequences: 100% identity, 100% coverage
- Spacer sequences: 90% identity, 90% coverage

Uses MMseqs2 for sequence clustering with proper handling of 
reverse complement sequences.
"""

from .cluster_crispr import (
    cluster_repeats,
    cluster_spacers,
    canonicalize_sequence,
    reverse_complement,
    check_mmseqs_installed,
)

__all__ = [
    'cluster_repeats',
    'cluster_spacers', 
    'canonicalize_sequence',
    'reverse_complement',
    'check_mmseqs_installed',
]
