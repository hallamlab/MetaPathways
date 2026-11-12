"""Ordered explorer categories; one accessible optgroup per branch."""
TABLE_GROUPS = {
    'Overview / samples and genomes': ['samples','entities','contig_mags'],
    'Features / contigs and genes': ['contigs','feature_explorer','orf_explorer','orf_groups'],
    'Features / RNA': ['rna_explorer','rna_taxonomy_explorer','orf_rna_explorer'],
    'Annotations / function and taxonomy': ['annotation_explorer','annotation_terms'],
    'MAGs / features': ['mag_orf_explorer','mag_rna_explorer','mag_gene_explorer'],
    'Pathways / inference and genes': ['pathway_explorer','pathway_gene_explorer'],
    'Abundance / measurements': ['abundance_explorer','rna_abundance_explorer'],
    'Provenance / files and diagnostics': ['sources','files','execution','issues'],
}
