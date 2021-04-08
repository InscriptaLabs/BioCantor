"""
Biotypes are types of genes and transcripts, as defined by NCBI (INSDC) and Sequence Ontology.
"""

from inscripta.biocantor.util.enum import HasMemberMixin


Biotype = HasMemberMixin(
    value="Biotype",
    names=[
        ["protein_coding", 0],
        ["protein-coding", 0],
        ["mRNA", 0],
        ["ncRNA", 1],
        ["misc_RNA", 2],
        ["miscRNA", 2],
        ["tRNA", 3],
        ["rRNA", 4],
        ["pseudogene", 5],
        ["pseudo", 5],
        ["lncRNA", 6],
        ["lnc_RNA", 6],
        ["snoRNA", 7],
        ["V_gene_segment", 8],
        ["C_gene_segment", 9],
        ["J_gene_segment", 10],
        ["D_gene_segment", 11],
        ["primary_transcript", 12],
        ["miRNA", 13],
        ["transcript", 14],
        ["SRP_RNA", 15],
        ["telomerase_RNA", 16],
        ["tmRNA", 17],
        ["RNase_MRP_RNA", 18],
        ["Y_RNA", 19],
        ["antisense_RNA", 20],
        ["scRNA", 21],
        ["snRNA", 22],
        ["vault_RNA", 23],
        ["J_segment", 24],
        ["C_region", 25],
        ["V_segment", 26],
        ["D_segment", 27],
        ["RNase_P_RNA", 28],
        ["transcribed_pseudogene", 29],
        ["ncRNA_pseudogene", 30],
        ["C_region_pseudogene", 31],
        ["J_segment_pseudogene", 32],
        ["V_segment_pseudogene", 33],
        ["guide_RNA", 34],
    ],
)


UNKNOWN_BIOTYPE = "unspecified"
