import json

import pytest
from Bio import SeqIO
from inscripta.biocantor.gene.biotype import Biotype
from inscripta.biocantor.gene.cds_frame import CDSFrame
from inscripta.biocantor.io.gff3.exc import GFF3FastaException
from inscripta.biocantor.io.exc import DuplicateSequenceException
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
                "Dbxref": ["Ensembl:ENST00000425520.2", "Genbank:NM_001080998.2", "GeneID:441581", "HGNC:HGNC:33518"],
                "gbkey": ["mRNA"],
                "gene": ["FRG2B"],
                "tag": ["MANE Select"],
            },
            "is_primary_tx": False,
            "transcript_id": "NM_001080998.2",
            "protein_id": "NP_001074467.1",
            "product": "protein FRG2-like-1",
            "transcript_symbol": None,
            "transcript_type": Biotype.protein_coding,
            "sequence_name": "NC_000010.11",
            "sequence_guid": None,
            "transcript_interval_guid": None,
            "transcript_guid": None,
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
        with open(test_data_dir / "SGCE.json") as fh:
            assert AnnotationCollectionModel.Schema().load(json.load(fh)) == rec.annotation

    def test_parse_peg10(self, test_data_dir):
        """
        PEG10 is a gene with a -1 frameshift in one isoform, that is parsed using the RefSeq parser.
        """
        gff = test_data_dir / "PEG10_minus1frameshift.gff3"
        rec = list(parse_standard_gff3(gff))[0]
        with open(test_data_dir / "PEG10_minus1frameshift.json") as fh:
            assert AnnotationCollectionModel.Schema().load(json.load(fh)) == rec.annotation

    def test_parse_insc1006(self, test_data_dir):
        """INSC1006_chrI is a 4-gene manually built file from INSC1006. It uses the default naive parser."""
        gff = test_data_dir / "INSC1006_chrI.gff3"
        recs = list(parse_standard_gff3(gff))
        c = recs[0].annotation
        with open(test_data_dir / "INSC1006_chrI.json") as fh:
            assert AnnotationCollectionModel.Schema().load(json.load(fh)) == c

    def test_parse_non_gene_locus_tag(self, test_data_dir):
        """Ran into bug in which non-gene feature with locus tag led to error, ensure handled correctly. Using cases
        from S288C gff3 that surfaced the issue."""
        gff = test_data_dir / "feature_test_non_gene_locus_tag.gff"
        recs = list(parse_standard_gff3(gff))
        features = list(recs[0].annotation.to_annotation_collection())
        assert features[0].to_dict()["locus_tag"] is None
        assert features[1].to_dict()["locus_tag"] == "YFL007W"

    @pytest.mark.parametrize(
        "gff3,json_file", [("feature_test_1.gff", "feature_test_1.json"), ("feature_test_2.gff", "feature_test_2.json")]
    )
    def test_parse_feature_tests(self, test_data_dir, gff3, json_file):
        recs = list(parse_standard_gff3(test_data_dir / gff3))
        c = recs[0].annotation
        with open(test_data_dir / json_file) as fh:
            assert AnnotationCollectionModel.Schema().load(json.load(fh)) == c

    def test_transcript_inference(self, test_data_dir):
        recs = list(parse_standard_gff3(test_data_dir / "feature_test_1.gff"))
        c = recs[0].annotation.to_annotation_collection()
        # 4 total genes
        assert len(c.genes) == 4
        # two different types of pseudogene transcripts were inferred, one without exons and one with exons
        # one gene with an invalid biotype who was set to None
        assert len([x for x in c.genes if x.gene_type == Biotype.pseudogene]) == 2
        assert len([x for x in c.genes if x.gene_type == Biotype.lncRNA]) == 1
        invalid_biotype = [x for x in c.genes if not x.gene_type]
        assert len(invalid_biotype) == 1
        assert list(invalid_biotype[0].qualifiers["provided_biotype"])[0] == "invalid"

    def test_direct_cds_exon(self, test_data_dir):
        recs = list(parse_standard_gff3(test_data_dir / "gene_cds_direct_child.gff3"))
        c = recs[0].annotation
        with open(test_data_dir / "gene_cds_direct_child.json") as fh:
            assert AnnotationCollectionModel.Schema().load(json.load(fh)) == c

    def test_cds_only(self, test_data_dir):
        """Some prokaryotic GFFs (e.g. prokka) do not have gene features, only CDS"""
        recs = list(parse_standard_gff3(test_data_dir / "test_cds_only.gff"))
        c = recs[0].annotation
        with open(test_data_dir / "test_cds_only.json") as fh:
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
        """INSC1003.gff3 has no sequences"""
        gff3_without_fasta = test_data_dir / "INSC1003.gff3"
        with pytest.raises(GFF3FastaException):
            _ = list(
                ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_gff3_embedded_fasta(gff3_without_fasta))
            )

    def test_duplicate_sequence(self, test_data_dir):
        fasta = test_data_dir / "INSC1003_extra_contig_duplicate.fa"
        gff3 = test_data_dir / "INSC1003.gff3"
        gff3_with_fasta = test_data_dir / "INSC1003_embedded_extra_contig_duplicate.gff3"
        with pytest.raises(DuplicateSequenceException):
            _ = list(parse_gff3_fasta(gff3, fasta))
        with pytest.raises(DuplicateSequenceException):
            _ = list(parse_gff3_embedded_fasta(gff3_with_fasta))


class TestGff3ToModel:
    """
    Test Prokaryotic genbanks. These have no mRNA feature, and they are inferred.

    The test case ``INSC1003.gff3`` is a subset of a Prokka annotation, with 7 total genes: 5 protein coding,
    1 tRNA and 1 rRNA. The tRNA and rRNA are not real, but the CDS are.
    """

    def test_gff3(self, test_data_dir):
        gff = test_data_dir / "INSC1003.gff3"
        fasta = test_data_dir / "INSC1003.fa"
        recs = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_gff3_fasta(gff, fasta)))
        for rec in recs:
            for gene in rec.genes:
                for tx in gene.transcripts:
                    # no frameshifts
                    if tx.is_coding:
                        assert "*" not in tx.get_protein_sequence()[:-1]

    def test_gff3_with_fa_extra_contig(self, test_data_dir):
        """Handle FASTA with sequences not seen in the GFF3"""
        gff3_without_fasta = test_data_dir / "INSC1003.gff3"
        fasta = test_data_dir / "INSC1003_extra_contig.fa"
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
        gff3 = test_data_dir / "INSC1003_embedded_extra_contig.gff3"
        recs = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_gff3_embedded_fasta(gff3)))
        assert len(recs) == 2
        assert recs[1].is_empty
        assert recs[0].sequence_name
        assert recs[0].sequence
        assert recs[1].sequence
        assert recs[1].sequence_name == "extraseq"
