"""
Biotypes are types of genes and transcripts, as defined by NCBI and Sequence Ontology.
"""

from inscripta.biocantor.util.enum import HasMemberMixin


Biotype = HasMemberMixin(
    value="Biotype",
    names=[
        ["protein_coding", 0],
        ["protein-coding", 1],
        ["ncRNA", 2],
        ["misc_RNA", 3],
        ["miscRNA", 3],
        ["tRNA", 4],
        ["rRNA", 5],
        ["pseudogene", 6],
        ["pseudo", 6],
        ["lncRNA", 7],
        ["lnc_RNA", 7],
        ["snoRNA", 8],
        ["V_gene_segment", 9],
        ["C_gene_segment", 10],
        ["J_gene_segment", 11],
        ["D_gene_segment", 12],
        ["primary_transcript", 13],
        ["miRNA", 14],
        ["transcript", 15],
        ["SRP_RNA", 16],
        ["telomerase_RNA", 17],
        ["tmRNA", 18],
        ["RNase_MRP_RNA", 19],
        ["Y_RNA", 20],
        ["antisense_RNA", 21],
        ["scRNA", 22],
        ["snRNA", 23],
        ["vault_RNA", 24],
        ["J_segment", 25],
        ["C_region", 26],
        ["V_segment", 27],
        ["D_segment", 28],
        ["RNase_P_RNA", 29],
        ["transcribed_pseudogene", 30],
    ],
)
