import json

import pytest
from Bio import SeqIO

from inscripta.biocantor.gene.biotype import Biotype
from inscripta.biocantor.gene.cds import CDSFrame
from inscripta.biocantor.io.gff3.exc import GFF3FastaException
from inscripta.biocantor.io.gff3.parser import (
    parse_gff3_embedded_fasta,
    parse_gff3_fasta,
    ParsedAnnotationRecord,
    parse_standard_gff3,
    extract_seqrecords_from_gff3_fasta,
)
from inscripta.biocantor.io.models import AnnotationCollectionModel
from inscripta.biocantor.location.strand import Strand


class TestGff3Parser:
    """Test GFF3 parsing"""

    def test_transitive(self, test_data_dir):
        """Test transitive loops"""
        for gff in ["SGCE.gff3", "FRG2B.gff3", "PEG10_minus1frameshift.gff3", "INSC1006_chrI.gff3"]:
            gff = test_data_dir / gff
            recs = list(parse_standard_gff3(gff))
            c = recs[0].to_annotation_collection()
            assert (
                AnnotationCollectionModel.Schema().load(c.to_dict()).to_annotation_collection().to_dict() == c.to_dict()
            )

    def test_parse_frg2b(self, test_data_dir):
        """
        FRG2B is an example of a standard coding gene.
        """
        gff = test_data_dir / "FRG2B.gff3"
        recs = list(parse_standard_gff3(gff))
        assert len(recs) == 1
        rec = recs[0]
        genes = rec.annotation.genes
        assert len(genes) == 1
        assert rec.annotation.feature_collections is None
        gene = genes[0]
        txs = gene.transcripts
        assert len(txs) == 1
        tx = txs[0]
        assert tx.__dict__ == {
            "exon_starts": [133623894, 133625921, 133626228, 133626564],
            "exon_ends": [133625604, 133625999, 133626303, 133626795],
            "strand": Strand.MINUS,
            "cds_starts": [133625098, 133625921, 133626228, 133626564],
            "cds_ends": [133625604, 133625999, 133626303, 133626742],
            "cds_frames": [CDSFrame.ONE, CDSFrame.ONE, CDSFrame.ONE, CDSFrame.ZERO],
            "qualifiers": {
                "Dbxref": ["Ensembl:ENST00000425520.2", "GeneID:441581", "Genbank:NM_001080998.2", "HGNC:HGNC:33518"],
                "gbkey": ["mRNA"],
                "gene": ["FRG2B"],
                "product": ["FSHD region gene 2 family member B"],
                "tag": ["MANE Select"],
                "transcript_id": ["NM_001080998.2"],
            },
            "is_primary_tx": False,
            "transcript_id": "NM_001080998.2",
            "protein_id": "NP_001074467.1",
            "transcript_symbol": None,
            "transcript_type": Biotype.protein_coding,
            "sequence_name": "NC_000010.11",
            "sequence_guid": None,
            "guid": None,
            "transcript_guid": None,
        }
        assert tx.exon_starts == [133623894, 133625921, 133626228, 133626564]
        assert tx.exon_ends == [133625604, 133625999, 133626303, 133626795]
        assert tx.cds_starts == [133625098, 133625921, 133626228, 133626564]
        assert tx.cds_ends == [133625604, 133625999, 133626303, 133626742]
        assert tx.qualifiers == {
            "Dbxref": ["Ensembl:ENST00000425520.2", "GeneID:441581", "Genbank:NM_001080998.2", "HGNC:HGNC:33518"],
            "gbkey": ["mRNA"],
            "gene": ["FRG2B"],
            "product": ["FSHD region gene 2 family member B"],
            "tag": ["MANE Select"],
            "transcript_id": ["NM_001080998.2"],
        }

    def test_parse_sgce(self, test_data_dir):
        """
        SGCE is an example of a protein coding gene with multiple isoforms, that is parsed using the RefSeq parser.
        """
        gff = test_data_dir / "SGCE.gff3"
        recs = list(parse_standard_gff3(gff))
        assert len(recs) == 1
        rec = recs[0]
        genes = rec.annotation.genes
        assert len(genes) == 1
        gene = genes[0]
        txs = gene.transcripts
        assert len(txs) == 22
        with open("tests/data/SGCE.json") as fh:
            assert AnnotationCollectionModel.Schema().load(json.load(fh)) == rec.annotation

    def test_parse_peg10(self, test_data_dir):
        """
        PEG10 is a gene with a -1 frameshift in one isoform, that is paresd using the RefSeq parser.
        """
        gff = test_data_dir / "PEG10_minus1frameshift.gff3"
        recs = list(parse_standard_gff3(gff))
        c = recs[0].to_annotation_collection()
        assert c.to_dict() == {
            "genes": [
                {
                    "transcripts": [
                        {
                            "exon_starts": (94656369, 94663333),
                            "exon_ends": (94656580, 94669695),
                            "strand": "PLUS",
                            "cds_starts": (94663556,),
                            "cds_ends": (94664534,),
                            "cds_frames": ["ZERO"],
                            "qualifiers": {
                                "Dbxref": ["GeneID:23089", "Genbank:NM_001040152.2", "HGNC:HGNC:14005", "MIM:609810"],
                                "gbkey": ["mRNA"],
                                "gene": ["PEG10"],
                                "product": ["paternally expressed 10, transcript variant 1"],
                                "transcript_id": ["NM_001040152.2"],
                            },
                            "is_primary_tx": False,
                            "transcript_id": "NM_001040152.2",
                            "transcript_symbol": None,
                            "transcript_type": "protein_coding",
                            "sequence_name": "chr7",
                            "sequence_guid": None,
                            "protein_id": "NP_001035242.1",
                            "transcript_guid": None,
                        },
                        {
                            "exon_starts": (94656324, 94663333),
                            "exon_ends": (94656591, 94669695),
                            "strand": "PLUS",
                            "cds_starts": (94656586, 94663333, 94664512),
                            "cds_ends": (94656591, 94664513, 94665682),
                            "cds_frames": ["ZERO", "TWO", "ZERO"],
                            "qualifiers": {
                                "Dbxref": ["GeneID:23089", "Genbank:NM_001172437.2", "HGNC:HGNC:14005", "MIM:609810"],
                                "gbkey": ["mRNA"],
                                "gene": ["PEG10"],
                                "product": ["paternally expressed 10, transcript variant 2"],
                                "transcript_id": ["NM_001172437.2"],
                            },
                            "is_primary_tx": False,
                            "transcript_id": "NM_001172437.2",
                            "transcript_symbol": None,
                            "transcript_type": "protein_coding",
                            "sequence_name": "chr7",
                            "sequence_guid": None,
                            "protein_id": "NP_001165908.1",
                            "transcript_guid": None,
                        },
                        {
                            "exon_starts": (94656369, 94663333),
                            "exon_ends": (94656591, 94669695),
                            "strand": "PLUS",
                            "cds_starts": (94656586, 94663333),
                            "cds_ends": (94656591, 94664534),
                            "cds_frames": ["ZERO", "TWO"],
                            "qualifiers": {
                                "Dbxref": ["GeneID:23089", "Genbank:NM_001172438.3", "HGNC:HGNC:14005", "MIM:609810"],
                                "gbkey": ["mRNA"],
                                "gene": ["PEG10"],
                                "product": ["paternally expressed 10, transcript variant 2"],
                                "transcript_id": ["NM_001172438.3"],
                            },
                            "is_primary_tx": False,
                            "transcript_id": "NM_001172438.3",
                            "transcript_symbol": None,
                            "transcript_type": "protein_coding",
                            "sequence_name": "chr7",
                            "sequence_guid": None,
                            "protein_id": "NP_001165909.1",
                            "transcript_guid": None,
                        },
                        {
                            "exon_starts": (94656324, 94663333),
                            "exon_ends": (94656580, 94669695),
                            "strand": "PLUS",
                            "cds_starts": (94663454, 94664512),
                            "cds_ends": (94664513, 94665682),
                            "cds_frames": ["ZERO", "ZERO"],
                            "qualifiers": {
                                "Dbxref": ["GeneID:23089", "Genbank:NM_001184961.1", "HGNC:HGNC:14005", "MIM:609810"],
                                "gbkey": ["mRNA"],
                                "gene": ["PEG10"],
                                "product": ["paternally expressed 10, transcript variant 1"],
                                "tag": ["RefSeq Select"],
                                "transcript_id": ["NM_001184961.1"],
                            },
                            "is_primary_tx": False,
                            "transcript_id": "NM_001184961.1",
                            "transcript_symbol": None,
                            "transcript_type": "protein_coding",
                            "sequence_name": "chr7",
                            "sequence_guid": None,
                            "protein_id": "NP_001171890.1",
                            "transcript_guid": None,
                        },
                        {
                            "exon_starts": (94656369, 94663333),
                            "exon_ends": (94656580, 94669695),
                            "strand": "PLUS",
                            "cds_starts": (94663454,),
                            "cds_ends": (94664534,),
                            "cds_frames": ["ZERO"],
                            "qualifiers": {
                                "Dbxref": ["GeneID:23089", "Genbank:NM_001184962.2", "HGNC:HGNC:14005", "MIM:609810"],
                                "gbkey": ["mRNA"],
                                "gene": ["PEG10"],
                                "product": ["paternally expressed 10, transcript variant 1"],
                                "transcript_id": ["NM_001184962.2"],
                            },
                            "is_primary_tx": False,
                            "transcript_id": "NM_001184962.2",
                            "transcript_symbol": None,
                            "transcript_type": "protein_coding",
                            "sequence_name": "chr7",
                            "sequence_guid": None,
                            "protein_id": "NP_001171891.1",
                            "transcript_guid": None,
                        },
                        {
                            "exon_starts": (94656324, 94663333),
                            "exon_ends": (94656580, 94669695),
                            "strand": "PLUS",
                            "cds_starts": (94663556, 94664512),
                            "cds_ends": (94664513, 94665682),
                            "cds_frames": ["ZERO", "ZERO"],
                            "qualifiers": {
                                "Dbxref": ["GeneID:23089", "Genbank:NM_015068.3", "HGNC:HGNC:14005", "MIM:609810"],
                                "gbkey": ["mRNA"],
                                "gene": ["PEG10"],
                                "product": ["paternally expressed 10, transcript variant 1"],
                                "transcript_id": ["NM_015068.3"],
                            },
                            "is_primary_tx": False,
                            "transcript_id": "NM_015068.3",
                            "transcript_symbol": None,
                            "transcript_type": "protein_coding",
                            "sequence_name": "chr7",
                            "sequence_guid": None,
                            "protein_id": "NP_055883.2",
                            "transcript_guid": None,
                        },
                    ],
                    "gene_id": None,
                    "gene_symbol": None,
                    "gene_type": "protein_coding",
                    "locus_tag": None,
                    "qualifiers": {
                        "Dbxref": ["GeneID:23089", "HGNC:HGNC:14005", "MIM:609810"],
                        "description": ["paternally expressed 10"],
                        "gbkey": ["Gene"],
                        "gene": ["PEG10"],
                        "gene_biotype": ["protein_coding"],
                        "gene_synonym": ["EDR", "HB-1", "Mar2", "Mart2", "MEF3L", "RGAG3", "RTL2", "SIRH1"],
                    },
                    "sequence_name": "chr7",
                    "sequence_guid": None,
                    "gene_guid": None,
                }
            ],
            "feature_collections": [],
            "name": None,
            "qualifiers": None,
            "sequence_name": "chr7",
            "sequence_guid": None,
            "start": 94656324,
            "end": 94669695,
            "completely_within": None,
        }

    def test_parse_insc1006(self, test_data_dir):
        """INSC1006_chrI is a 4-gene manually built file from INSC1006. It uses the default naive parser."""
        gff = test_data_dir / "INSC1006_chrI.gff3"
        recs = list(parse_standard_gff3(gff))
        c = recs[0].annotation
        with open("tests/data/INSC1006_chrI.json") as fh:
            assert AnnotationCollectionModel.Schema().load(json.load(fh)) == c


class TestGff3FastaParser:
    """Test GFF3 + FASTA parsing"""

    def test_gff3_with_fa(self, test_data_dir):
        """Parse the FASTA and compare to what we get from GenBank"""
        gff3_with_fa = test_data_dir / "INSC1006_chrI.gff3"
        gbk = str(test_data_dir / "INSC1006_chrI.gbff")
        with open(gff3_with_fa, "r") as fh:
            recs = extract_seqrecords_from_gff3_fasta(fh)
        gbk_recs = list(SeqIO.parse(gbk, format="genbank"))
        for gff3_rec, gbk_rec in zip(recs, gbk_recs):
            assert gff3_rec.seq == gbk_rec.seq
            assert gff3_rec.id == gbk_rec.id
            # BioPython GenBank parser checks for ACCESSION.VESRION, FASTA parser does not.
            assert gff3_rec.name.split(".")[0] == gbk_rec.name

    def test_no_fa_exception(self, test_data_dir):
        """Inscripta_BL21.gff3 has no sequences"""
        gff3_without_fasta = test_data_dir / "Inscripta_BL21.gff3"
        with pytest.raises(GFF3FastaException):
            _ = list(
                ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_gff3_embedded_fasta(gff3_without_fasta))
            )


class TestGff3ToModel:
    """
    Test Prokaryotic genbanks. These have no mRNA feature, and they are inferred.

    The test case ``Inscripta_BL21.gff3`` is a subset of a Prokka annotation, with 7 total genes: 5 protein coding,
    1 tRNA and 1 rRNA. The tRNA and rRNA are not real, but the CDS are.
    """

    def test_gff3(self, test_data_dir):
        gff = test_data_dir / "Inscripta_BL21.gff3"
        fasta = test_data_dir / "Inscripta_BL21.fa"
        recs = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_gff3_fasta(gff, fasta)))
        for rec in recs:
            for gene in rec.genes:
                for tx in gene.transcripts:
                    # no frameshifts
                    if tx.is_coding:
                        assert "*" not in tx.get_protein_sequence()[:-1]

    def test_gff3_with_fa_extra_contig(self, test_data_dir):
        """Handle FASTA with sequences not seen in the GFF3"""
        gff3_without_fasta = test_data_dir / "Inscripta_BL21.gff3"
        fasta = test_data_dir / "Inscripta_BL21_extra_contig.fa"
        recs = list(
            ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_gff3_fasta(gff3_without_fasta, fasta))
        )
        assert len(recs) == 2
        assert recs[1].is_empty
        assert recs[0].sequence_name
        assert recs[0].sequence
        assert recs[1].sequence
        assert recs[1].sequence_name == "extraseq"

    def test_gff3_with_embedded_fa_extra_contig(self, test_data_dir):
        """Handle FASTA with sequences not seen in the GFF3"""
        gff3 = test_data_dir / "Inscripta_BL21_embedded_extra_contig.gff3"
        recs = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_gff3_embedded_fasta(gff3)))
        assert len(recs) == 2
        assert recs[1].is_empty
        assert recs[0].sequence_name
        assert recs[0].sequence
        assert recs[1].sequence
        assert recs[1].sequence_name == "extraseq"
