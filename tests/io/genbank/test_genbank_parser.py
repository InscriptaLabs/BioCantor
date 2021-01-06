from uuid import UUID

import pytest
from inscripta.biocantor.gene.biotype import Biotype
from inscripta.biocantor.gene.cds import CDSFrame
from inscripta.biocantor.io.genbank.exc import GenBankValidationError
from inscripta.biocantor.io.genbank.parser import parse_genbank, GenBankParserType
from inscripta.biocantor.io.parser import ParsedAnnotationRecord
from inscripta.biocantor.location.location_impl import SingleInterval, CompoundInterval, Strand
from inscripta.biocantor.util.hashing import digest_object


class TestEukaryoticGenbankParser:
    """
    The test case ``INSC1006_chrI.gbff`` has 5 genes. Each tests something different:

    GI526_G0000001: A non-coding gene.
    GI526_G0000002: A eukaryotic style coding gene. Has 5' and 3' UTRs.
    GI526_G0000003: A prokaryotic style coding gene. Has no mRNA feature and must be inferred.
    GI526_G0000004: A coding gene with no CDS feature; gets parsed as a non-coding gene.
    GI526_G0000005: Gene feature only; should not be parsed.

    """

    gbk = "INSC1006_chrI.gbff"

    def test_parse_genbank_metadata(self, test_data_dir):
        gbk = test_data_dir / self.gbk
        with open(gbk, "r") as fh:
            parsed = list(parse_genbank(fh))[0]

        assert parsed.annotation.feature_collections is None
        assert parsed.annotation.sequence_name == "CM021111.1"
        assert len(parsed.annotation.genes) == 4
        assert {x.gene_symbol for x in parsed.annotation.genes if x.gene_symbol} == {"GDH3", "BDH2", "BDH1"}
        assert {x.locus_tag for x in parsed.annotation.genes} == {
            "GI526_G0000001",
            "GI526_G0000002",
            "GI526_G0000003",
            "GI526_G0000004",
        }


class TestSplicedGenbankParser:
    """Test Spliced GenBank parsing"""

    gbk = "test_spliced.gbff"

    def test_parse_genbank_metadata(self, test_data_dir):
        gbk = test_data_dir / self.gbk
        with open(gbk, "r") as fh:
            parsed = list(parse_genbank(fh))[0]

        assert {x.gene_symbol for x in parsed.annotation.genes if x.gene_symbol} == {"MPT5"}


class TestProkaryoticGenBankParser:
    """
    Test Prokaryotic GenBanks. These have no mRNA feature, and they are inferred.

    The test case ``Inscripta_BL21.gbk`` is a subset of a Prokka annotation, with 7 total genes: 5 protein coding,
    1 tRNA and 1 rRNA. The tRNA and rRNA are not real, but the CDS are.
    """

    gbk = "Inscripta_BL21.gbk"

    def test_parse_genbank_metadata(self, test_data_dir):
        gbk = test_data_dir / self.gbk
        with open(gbk, "r") as fh:
            parsed = list(parse_genbank(fh))[0]

        assert {x.gene_symbol for x in parsed.annotation.genes if x.gene_symbol} == {
            "thrA",
            "thrB",
            "thrC",
            "yaaA",
            "yaaX",
        }
        assert {x.locus_tag for x in parsed.annotation.genes} == {
            "FEPOIHMA_00001",
            "FEPOIHMA_00002",
            "FEPOIHMA_00003",
            "FEPOIHMA_00004",
            "FEPOIHMA_00005",
            "FEPOIHMA_02080",
            "FEPOIHMA_02457",
        }
        assert [x.gene_type for x in parsed.annotation.genes] == [Biotype.protein_coding] * 5 + [
            Biotype.tRNA,
            Biotype.rRNA,
        ]


class TestGenbank:
    """
    The test case ``INSC1006_chrI.gbff`` has 5 genes. Each tests something different:

    GI526_G0000001: A non-coding gene.
    GI526_G0000002: A eukaryotic style coding gene. Has 5' and 3' UTRs.
    GI526_G0000003: A prokaryotic style coding gene. Has no mRNA feature and must be inferred.
    GI526_G0000004: A coding gene with no CDS feature; gets parsed as a non-coding gene.
    GI526_G0000005: Gene feature only; should not be parsed.

    """

    gbk = "INSC1006_chrI.gbff"

    def test_genbank_to_model(self, test_data_dir):
        gbk = test_data_dir / self.gbk
        with open(gbk, "r") as fh:
            parsed = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(fh)))[0]

        assert len(parsed.genes) == 4
        assert all(
            gene.transcripts[0].location.parent.id == gene.sequence_name == "CM021111.1" for gene in parsed.genes
        )
        assert digest_object(parsed) == UUID("6845b856-6135-76c5-8f9c-b127bddb2b9f")
        assert not parsed.genes[0].transcripts[0].is_coding
        assert parsed.genes[1].transcripts[0].is_coding
        assert parsed.genes[2].transcripts[0].is_coding
        assert not parsed.genes[3].transcripts[0].is_coding
        # has UTR
        assert parsed.genes[1].transcripts[0].location != parsed.genes[1].transcripts[0].cds.location
        # does not have UTR
        assert parsed.genes[2].transcripts[0].location == parsed.genes[2].transcripts[0].cds.location

        # validate positions; gene always has + strand location
        for gene, expected_gene_loc, expected_tx_loc in zip(
            parsed,
            ["16174-18079:+", "37461-39103:+", "39518-40772:+", "41085-42503:+"],
            ["16174-18079:-", "37461-39103:+", "39518-40772:+", "41085-42503:+"],
        ):
            assert str(gene.location) == expected_gene_loc
            assert str(gene.transcripts[0].location) == expected_tx_loc

        # has UTR
        assert parsed.genes[1].transcripts[0].cds_size == 1374
        assert len(parsed.genes[1].transcripts[0]) == 1642
        # does not have UTR
        assert parsed.genes[2].transcripts[0].cds_size == 1254
        assert len(parsed.genes[2].transcripts[0]) == 1254

    @pytest.mark.parametrize(
        "seq_id,start,end,coding_only,completely_within,locus_tags",
        [
            (
                # empty filter; pass through
                None,
                None,
                None,
                False,
                True,
                {"GI526_G0000001", "GI526_G0000002", "GI526_G0000003", "GI526_G0000004"},
            ),
            (
                # just Seq ID is pass through
                "CM021111.1",
                None,
                None,
                False,
                True,
                {"GI526_G0000001", "GI526_G0000002", "GI526_G0000003", "GI526_G0000004"},
            ),
            (
                # restrict by strict range
                "CM021111.1",
                10000,
                37600,
                False,
                True,
                {"GI526_G0000001"},
            ),
            (
                # restrict by not strict range
                "CM021111.1",
                10000,
                37600,
                False,
                False,
                {"GI526_G0000001", "GI526_G0000002"},
            ),
            (
                # restrict by coding
                None,
                None,
                None,
                True,
                True,
                {"GI526_G0000002", "GI526_G0000003"},
            ),
            (
                # exact match
                "CM021111.1",
                16174,
                18079,
                False,
                True,
                {"GI526_G0000001"},
            ),
        ],
    )
    def test_interval_filtering(self, test_data_dir, seq_id, start, end, coding_only, completely_within, locus_tags):
        gbk = test_data_dir / self.gbk
        with open(gbk, "r") as fh:
            annot_collection = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(fh)))[0]

        result = annot_collection.query_by_position(
            start=start,
            end=end,
            coding_only=coding_only,
            completely_within=completely_within,
        )
        assert all(x.locus_tag in locus_tags for x in result.genes)

    @pytest.mark.parametrize(
        "seq_id,start,end,coding_only,completely_within,locus_tags",
        [
            (
                # restrict by range and coding
                "CM021111.1",
                10000,
                37600,
                True,
                True,
                None,
            ),
            (
                # restrict to exactly 1bp off
                "CM021111.1",
                10000,
                18078,
                True,
                True,
                None,
            ),
            (
                # restrict to exactly 1bp off
                "CM021111.1",
                16175,
                18079,
                True,
                True,
                None,
            ),
        ],
    )
    def test_empty_interval_filtering(
        self, test_data_dir, seq_id, start, end, coding_only, completely_within, locus_tags
    ):
        gbk = test_data_dir / self.gbk
        with open(gbk, "r") as fh:
            annot_collection = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(fh)))[0]

        result = annot_collection.query_by_position(
            start=start,
            end=end,
            coding_only=coding_only,
            completely_within=completely_within,
        )
        assert result.is_empty

    @pytest.mark.parametrize(
        "identifiers,locus_tags",
        [
            (
                # locus tag
                ["GI526_G0000001", "GI526_G0000002", "GI526_G0000003", "GI526_G0000004"],
                {"GI526_G0000001", "GI526_G0000002", "GI526_G0000003", "GI526_G0000004"},
            ),
            (
                # locus tag
                ["GI526_G0000001", "GI526_G0000002", "GI526_G0000003"],
                {"GI526_G0000001", "GI526_G0000002", "GI526_G0000003"},
            ),
            (
                # gene ID
                ["GDH3"],
                {"GI526_G0000002"},
            ),
        ],
    )
    def test_identifier_filtering(self, test_data_dir, identifiers, locus_tags):
        gbk = test_data_dir / self.gbk
        with open(gbk, "r") as fh:
            annot_collection = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(fh)))[0]

        result = annot_collection.query_by_feature_identifier(identifiers)
        assert all(x.locus_tag in locus_tags for x in result.genes)

    @pytest.mark.parametrize(
        "seq_id,identifiers",
        [("CM021111.1", ["GI526_G1000001"])],
    )
    def test_identifier_empty_filtering(self, test_data_dir, seq_id, identifiers):
        gbk = test_data_dir / self.gbk
        with open(gbk, "r") as fh:
            annot_collection = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(fh)))[0]

        result = annot_collection.query_by_feature_identifier(identifiers)
        assert result.is_empty


class TestSplicedGenbank:
    """Test Spliced GenBank Parsing"""

    gbk = "test_spliced.gbff"

    def test_parsing(self, test_data_dir):
        gbk = test_data_dir / self.gbk
        with open(gbk, "r") as fh:
            annot_collection = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(fh)))[0]

        spliced = annot_collection.genes[0]
        assert spliced.gene_symbol == "MPT5"
        assert spliced.location.reset_parent(None) == SingleInterval(1000, 4220, Strand.PLUS)
        assert spliced.get_primary_transcript().location.reset_parent(None) == CompoundInterval(
            [1000, 1643], [1003, 4220], Strand.PLUS
        )
        assert spliced.get_primary_cds().location.reset_parent(None) == CompoundInterval(
            [1000, 1643], [1003, 4220], Strand.PLUS
        )
        assert (
            str(spliced.get_primary_protein())
            == "MINNEPFPSADSASILTTSTSNNSLMSYNHQPQLSINSVQSLLEPVTPPPLGQMNNKRNHQKAHSLDLSGFNQFISSTQSPLALMNNTSTSNSANSFSPNP"
            "NAASNSTGLSASMANPPAILPLINEFDLEMDGPRRKSSHDFTVVAPSNSGVNTSSLIMETPSSSVTPAASLRNFSNSNNAASKCGVDNSSFGLSSSTSSSM"
            "VEISALPLRDLDYIKLATDQFGCRFLQKKLETPSESNMVRDLMYEQIKPFFLDLILDPFGNYLVQKLCDYLTAEQKTLLIQTIYPNVFQISINQYGTRSLQ"
            "KIIDTVDNEVQIDLIIKGFSQEFTSIEQVVTLINDLNGNHVIQKCIFKFSPSKFGFIIDAIVEQNNIITISTHKHGCCVLQKLLSVCTLQQIFKISVKIVQ"
            "FLPGLINDQFGNYIIQFLLDIKELDFYLLAELFNRLSNELCQLSCLKFSSNVVEKFIKKLFRIITGFIVNNNGGASQRTAVASDDVINASMNILLTTIDIF"
            "TVNLNVLIRDNFGNYALQTLLDVKNYSPLLAYNKNSNAIGQNSSSTLNYGNFCNDFSLKIGNLIVLTKELLPSIKTTSYAKKIKLKVKAYAEATGIPFTDI"
            "SPQVTAMSHNNLQTINNENKNPHNKNSHNHNHNHNHNHAHNNNNNNNQKSHTRHFSLPANAYHRRSNSSVTNNFSNQYAQDQKIHSPQQIMNFNQNAYPSM"
            "GAPSFNSQTNPPLVSHNSLQNFDNRQFANLMAHPNSAAPIHSFSSSNITNVNPNVSRGFKQPGFMMNETDKINANHFSPYSNANSQNFNESFVPRMQYQTE"
            "GANWDSSLSMKSQHIGQGPYNQVNMSRNASISNMPAMNTARTSDELQFTLP*"
        )


class TestProkaryoticGenBank:
    """
    Test Prokaryotic genbanks. These have no mRNA feature, and they are inferred.

    The test case ``Inscripta_BL21.gbk`` is a subset of a Prokka annotation, with 7 total genes: 5 protein coding,
    1 tRNA and 1 rRNA. The tRNA and rRNA are not real, but the CDS are.
    """

    gbk = "Inscripta_BL21.gbk"

    def test_parsing(self, test_data_dir):
        gbk = test_data_dir / self.gbk
        with open(gbk, "r") as fh:
            annot_collection = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(fh)))[0]

        for gene in annot_collection.genes:
            for tx in gene.transcripts:
                if tx.is_coding:
                    # check that protein equals what Prokka thought it was
                    assert {str(tx.get_protein_sequence())[:-1]} == tx.qualifiers["translation"]


class TestFrameGenBank:
    """
    Frame has to be inferred from GenBank files because there is no Frame/Phase field. There is a /codon_start field,
    which works for offset start codons, but cannot model indels or programmed frameshifts. This set of tests
    makes sure that GenBank parsing properly handles building the frames, especially on the negative strand.
    """

    def test_negative_strand(self, test_data_dir):
        """
        This file has the same transcript in three versions:
        1. No /codon_start (infer /codon_start=1)
        2. Explicit /codon_start (infer /codon_start=1)
        3. /codon_start=2. This leads to the first codon being removed.

        This transcript is extra funky because it has exactly 1 base of the start codon present on the 2nd, 1bp, exon.

        In this case, because the /codon_start=2 truncates the 2nd exon entirely, the frame ends up undisturbed
        in the main exon.
        """
        gbk = test_data_dir / "negative_strand_frame.gbk"
        with open(gbk, "r") as fh:
            annot_collection = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(fh)))[0]

        prot_str = (
            "MAHFKEYQVIGRRLPTESVPEPKLFRMRIFASNEVIAKSRYWYFLQKLHKVKKASGEIVSINQINEAHPTKVKNFGVWVRYDSRSGTHNMYKEIRD"
            "VSRVAAVETLYQDMAARHRARFRSIHILKVAEIEKTADVKRQYVKQFLTKDLKFPLPHRVQKSTKTFSYKRPSTFY*"
        )

        assert (
            str(annot_collection.genes[0].get_primary_protein())
            == str(annot_collection.genes[1].get_primary_protein())
            == prot_str
        )
        assert (
            annot_collection.genes[0].get_primary_cds().frames
            == annot_collection.genes[1].get_primary_cds().frames
            == [CDSFrame.ONE, CDSFrame.ZERO]
        )

        assert str(annot_collection.genes[2].get_primary_protein()) == prot_str[1:]
        assert annot_collection.genes[2].get_primary_cds().frames == [CDSFrame.TWO, CDSFrame.ONE]

    def test_positive_strand(self, test_data_dir):
        """
        This file has the same transcript in three versions:
        1. No /codon_start (infer /codon_start=1)
        2. Explicit /codon_start (infer /codon_start=1)
        3. /codon_start=2. This will lead to a different frame set and the first codon being removed.
        """
        gbk = test_data_dir / "positive_strand_frame.gbk"
        with open(gbk, "r") as fh:
            annot_collection = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(fh)))[0]

        prot_str = (
            "MSTEKILTPESQLKKTKAQQKTAEQIAAERAARKAANKEKRAIILERNAAYQKEYETAERNIIQAKRDAKAAGSYYVEAQHKLVFVVRIKGINKIPPKPRKV"
            "LQLLRLTRINSGTFVKVTKATLELLKLIEPYVAYGYPSYSTIRQLVYKRGFGKINKQRVPLSDNAIIEANLGKYGILSIDDLIHEIITVGPHFKQANNFLWP"
            "FKLSNPSGGWGVPRKFKHFIQGGSFGNREEFINKLVKAMN*"
        )

        assert (
            str(annot_collection.genes[0].get_primary_protein())
            == str(annot_collection.genes[1].get_primary_protein())
            == prot_str
        )
        assert (
            annot_collection.genes[0].get_primary_cds().frames
            == annot_collection.genes[1].get_primary_cds().frames
            == [CDSFrame.ZERO, CDSFrame.TWO, CDSFrame.ZERO]
        )

        assert (
            str(annot_collection.genes[2].get_primary_protein())
            == "CPLKNLDS*ISIEED*SSTKDCRTNCCRESCP*SR*QGKKSYYFGKKRRLPKGIRNC*KKHHSS*A*CQGCWFLLRRSSTQVGLRCQNQGY*QDST*AKKG"
            "STIAKIDKNQLWYIRQSYQGYFGTIEVD*TIRCLRLPILLYY*TIGLQERFR*DQQAKSSIVRQCYHRSQLG*VWYLVH*RFDSRNHHCWSTLQAS*QLFV"
            "AIQVVQPIWWLGCPKKVQAFHPRWFFR*P*RIHQ*IG*GYEL"
        )

        assert (
            str(annot_collection.genes[3].get_primary_protein())
            == "VH*KILTPESQLKKTKAQQKTAEQIAAERAARKA*QGKKSYYFGKKRRLPKGIRNC*KKHHSS*A*CQGCWFLLRRSSTQVGLRCQNQGY*QDST*AKKG"
            "STIAKIDKNQLWYIRQSYQGYFGTIEVD*TIRCLRLPILLYY*TIGLQERFR*DQQAKSSIVRQCYHRSQLG*VWYLVH*RFDSRNHHCWSTLQAS*QLF"
            "VAIQVVQPIWWLGCPKKVQAFHPRWFFR*P*RIHQ*IG*GYEL"
        )


class TestGenBankErrors:
    def test_locus_tag_unique(self, test_data_dir):
        """If using the default locus_tag way of detecting groupings, locus tags must be unique"""
        gbk = test_data_dir / "locus_tag_collision.gbk"
        with pytest.raises(GenBankValidationError):
            with open(gbk, "r") as fh:
                _ = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(fh)))[0]
        # works fine in sorted mode
        with open(gbk, "r") as fh:
            rec = list(
                ParsedAnnotationRecord.parsed_annotation_records_to_model(
                    parse_genbank(fh, gbk_type=GenBankParserType.SORTED)
                )
            )[0]
            assert rec.genes[0].locus_tag == rec.genes[1].locus_tag
