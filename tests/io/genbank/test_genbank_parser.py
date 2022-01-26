import json
from collections import OrderedDict

import pytest
from inscripta.biocantor.gene.biotype import Biotype
from inscripta.biocantor.gene.cds_frame import CDSFrame
from inscripta.biocantor.io.exc import StrandViolationWarning, DuplicateSequenceException
from inscripta.biocantor.io.genbank.exc import (
    GenBankLocusTagError,
    GenBankLocationException,
)
from inscripta.biocantor.io.genbank.parser import parse_genbank, GenBankParserType
from inscripta.biocantor.io.models import AnnotationCollectionModel
from inscripta.biocantor.io.parser import ParsedAnnotationRecord
from inscripta.biocantor.location.location_impl import SingleInterval, CompoundInterval, Strand


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

    The test case ``INSC1003.gbk`` is a subset of a Prokka annotation, with 7 total genes: 5 protein coding,
    1 tRNA and 1 rRNA. The tRNA and rRNA are not real, but the CDS are.
    """

    gbk = "INSC1003.gbk"

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
            gene.transcripts[0]._location.parent.id == gene.sequence_name == "CM021111.1" for gene in parsed.genes
        )
        assert not parsed.genes[0].transcripts[0].is_coding
        assert parsed.genes[1].transcripts[0].is_coding
        assert parsed.genes[2].transcripts[0].is_coding
        assert not parsed.genes[3].transcripts[0].is_coding
        # has UTR
        assert parsed.genes[1].transcripts[0]._location != parsed.genes[1].transcripts[0].cds._location
        # does not have UTR
        assert parsed.genes[2].transcripts[0]._location == parsed.genes[2].transcripts[0].cds._location

        # validate positions; gene always has + strand location
        for gene, expected_gene_loc, expected_tx_loc in zip(
            parsed,
            ["16174-18079:+", "37461-39103:+", "39518-40772:+", "41085-42503:+"],
            ["16174-18079:-", "37461-39103:+", "39518-40772:+", "41085-42503:+"],
        ):
            assert str(gene._location) == expected_gene_loc
            assert str(gene.transcripts[0]._location) == expected_tx_loc

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

        result = annot_collection.query_by_feature_identifiers(identifiers)
        assert all(x.locus_tag in locus_tags for x in result.genes)

    @pytest.mark.parametrize(
        "seq_id,identifiers",
        [("CM021111.1", ["GI526_G1000001"])],
    )
    def test_identifier_empty_filtering(self, test_data_dir, seq_id, identifiers):
        gbk = test_data_dir / self.gbk
        with open(gbk, "r") as fh:
            annot_collection = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(fh)))[0]

        result = annot_collection.query_by_feature_identifiers(identifiers)
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
        assert spliced._location.reset_parent(None) == SingleInterval(1000, 4220, Strand.PLUS)
        assert spliced.get_primary_transcript()._location.reset_parent(None) == CompoundInterval(
            [1000, 1643], [1003, 4220], Strand.PLUS
        )
        assert spliced.get_primary_cds()._location.reset_parent(None) == CompoundInterval(
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

    The test case ``INSC1003.gbk`` is a subset of a Prokka annotation, with 7 total genes: 5 protein coding,
    1 tRNA and 1 rRNA. The tRNA and rRNA are not real, but the CDS are.
    """

    gbk = "INSC1003.gbk"

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

        """
        gbk = test_data_dir / "negative_strand_frame.gbk"
        with open(gbk, "r") as fh:
            annot_collection = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(fh)))[0]

        prot_str = (
            "MAHFKEYQVIGRRLPTESVPEPKLFRMRIFASNEVIAKSRYWYFLQKLHKVKKASGEIVSINQINEAHPTKVKNFGVWVRYDSRSGTHNMYKEIRD"
            "VSRVAAVETLYQDMAARHRARFRSIHILKVAEIEKTADVKRQYVKQFLTKDLKFPLPHRVQKSTKTFSYKRPSTFY*"
        )

        shifted_prot_str = (
            "WLTLKNTKLLAVVCQLNLFQNQSCSE*ESLLQMKLLPSLVTGISCKSCTRLRRLLVKLFPSTKSTKLIQPRSRTSVSGLDTTPDLVL"
            "TICTRKSETSPELLPSKPYTKTWLPDTELDLDLFTS*RLLKLKRLLTSRDNTLSNF*PRT*NSHCLTESKNPPRLSPTRDLPLST"
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

        assert str(annot_collection.genes[2].get_primary_protein()) == shifted_prot_str
        assert annot_collection.genes[2].get_primary_cds().frames == [CDSFrame.ZERO, CDSFrame.ONE]

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
            == "CPLKKS*LLNLN*RRLKLNKRLQNKLLQRELPVKPLTRKKELLFWKETPLTKRNTKLLKETSFKLSVMPRLLVPTTSKLNTSWSSLSESRVLTRFHLSQER"
            "FYNC*D*QESTLVHSSKLPRLLWNY*S*LNHTLLTVTHPTLLLDNWSTREVSVRSTSKEFHCPTMLSSKPTWVSMVSCPLTI*FTKSSLLVHTSSKLTTFC"
            "GHSSCPTHLVVGVSQESSSISSKVVLSVTVKNSSINWLRL*T"
        )

        assert (
            str(annot_collection.genes[3].get_primary_protein())
            == "VH*KNLDS*ISIEED*SSTKDCRTNCCRESCP*SR*QGKKSYYFGKKRRLPKGIRNC*KKHHSS*A*CQGCWFLLRRSSTQVGLRCQNQGY*QDST*"
            "AKKGSTIAKIDKNQLWYIRQSYQGYFGTIEVD*TIRCLRLPILLYY*TIGLQERFR*DQQAKSSIVRQCYHRSQLG*VWYLVH*RFDSRNHHCWSTLQ"
            "AS*QLFVAIQVVQPIWWLGCPKKVQAFHPRWFFR*P*RIHQ*IG*GYEL"
        )


class TestGenBankErrors:
    def test_locus_tag_unique(self, test_data_dir):
        """If using the default locus_tag way of detecting groupings, locus tags must be unique"""
        gbk = test_data_dir / "locus_tag_collision.gbk"
        with pytest.raises(GenBankLocusTagError):
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


class TestGenBankFeatures:
    def test_parse_feature_test_2(self, test_data_dir):
        genbank = "feature_test_2.gbk"
        json_file = "feature_test_2_gbk.json"
        recs = list(parse_genbank(test_data_dir / genbank))
        c = recs[0].annotation
        assert len(c.feature_collections) == 2
        assert len(c.feature_collections[0].feature_intervals) == 1
        assert len(c.feature_collections[1].feature_intervals) == 3
        assert len(c.genes) == 1

        with open(test_data_dir / json_file) as fh:
            assert AnnotationCollectionModel.Schema().load(json.load(fh)) == c

    def test_parse_feature_test_3(self, test_data_dir):
        genbank = "feature_test_3.gbk"
        recs = list(parse_genbank(test_data_dir / genbank))
        c = recs[0].annotation
        assert len(c.feature_collections) == 4
        assert len(c.genes) == 1
        assert AnnotationCollectionModel.Schema().dump(c) == OrderedDict(
            [
                (
                    "feature_collections",
                    [
                        OrderedDict(
                            [
                                (
                                    "feature_intervals",
                                    [
                                        OrderedDict(
                                            [
                                                ("interval_starts", [17999]),
                                                ("interval_ends", [19000]),
                                                ("strand", "PLUS"),
                                                (
                                                    "qualifiers",
                                                    {
                                                        "name": ["overlapping_ncrna_opposite_strand"],
                                                        "note": ["not a joined interval"],
                                                    },
                                                ),
                                                ("sequence_name", "CM021111.1"),
                                                ("sequence_guid", None),
                                                ("feature_interval_guid", "dad6ee93-38c2-bf22-5111-d7c0fc0a5863"),
                                                ("feature_guid", None),
                                                ("feature_types", ["feature"]),
                                                ("feature_name", "overlapping_ncrna_opposite_strand"),
                                                ("feature_id", None),
                                                ("is_primary_feature", False),
                                            ]
                                        )
                                    ],
                                ),
                                ("feature_collection_name", "overlapping_ncrna_opposite_strand"),
                                ("feature_collection_id", None),
                                ("locus_tag", None),
                                ("feature_collection_type", None),
                                ("sequence_name", "CM021111.1"),
                                ("sequence_guid", None),
                                ("feature_collection_guid", "cc74d80d-d815-7b93-f105-8801c81cb164"),
                                (
                                    "qualifiers",
                                    {"name": ["overlapping_ncrna_opposite_strand"], "note": ["not a joined interval"]},
                                ),
                            ]
                        ),
                        OrderedDict(
                            [
                                (
                                    "feature_intervals",
                                    [
                                        OrderedDict(
                                            [
                                                ("interval_starts", [20000, 20500]),
                                                ("interval_ends", [20250, 20600]),
                                                ("strand", "PLUS"),
                                                ("qualifiers", {"name": ["joined_feature_plus_strand"]}),
                                                ("sequence_name", "CM021111.1"),
                                                ("sequence_guid", None),
                                                ("feature_interval_guid", "b048fa8c-e5e9-b903-42f4-ae7275a79506"),
                                                ("feature_guid", None),
                                                ("feature_types", ["feature"]),
                                                ("feature_name", "joined_feature_plus_strand"),
                                                ("feature_id", None),
                                                ("is_primary_feature", False),
                                            ]
                                        )
                                    ],
                                ),
                                ("feature_collection_name", "joined_feature_plus_strand"),
                                ("feature_collection_id", None),
                                ("locus_tag", None),
                                ("feature_collection_type", None),
                                ("sequence_name", "CM021111.1"),
                                ("sequence_guid", None),
                                ("feature_collection_guid", "e986d20b-2e0e-8e11-76d1-bd49796ef25c"),
                                ("qualifiers", {"name": ["joined_feature_plus_strand"]}),
                            ]
                        ),
                        OrderedDict(
                            [
                                (
                                    "feature_intervals",
                                    [
                                        OrderedDict(
                                            [
                                                ("interval_starts", [20550, 20749]),
                                                ("interval_ends", [20700, 21000]),
                                                ("strand", "MINUS"),
                                                (
                                                    "qualifiers",
                                                    {
                                                        "name": ["joined_feature_minus_strand"],
                                                        "note": ["overlaps last feature"],
                                                    },
                                                ),
                                                ("sequence_name", "CM021111.1"),
                                                ("sequence_guid", None),
                                                ("feature_interval_guid", "da172494-6af0-688e-5ff2-9bc583fec47d"),
                                                ("feature_guid", None),
                                                ("feature_types", ["feature"]),
                                                ("feature_name", "joined_feature_minus_strand"),
                                                ("feature_id", None),
                                                ("is_primary_feature", False),
                                            ]
                                        )
                                    ],
                                ),
                                ("feature_collection_name", "joined_feature_minus_strand"),
                                ("feature_collection_id", None),
                                ("locus_tag", None),
                                ("feature_collection_type", None),
                                ("sequence_name", "CM021111.1"),
                                ("sequence_guid", None),
                                ("feature_collection_guid", "46d4dd7a-6c99-0264-4def-6180624f84c5"),
                                (
                                    "qualifiers",
                                    {"name": ["joined_feature_minus_strand"], "note": ["overlaps last feature"]},
                                ),
                            ]
                        ),
                        OrderedDict(
                            [
                                (
                                    "feature_intervals",
                                    [
                                        OrderedDict(
                                            [
                                                ("interval_starts", [24999]),
                                                ("interval_ends", [26000]),
                                                ("strand", "MINUS"),
                                                ("qualifiers", {"name": ["unjoined_minus_strand"]}),
                                                ("sequence_name", "CM021111.1"),
                                                ("sequence_guid", None),
                                                ("feature_interval_guid", "a853506f-377f-8093-6f01-342c2c61be69"),
                                                ("feature_guid", None),
                                                ("feature_types", ["feature"]),
                                                ("feature_name", "unjoined_minus_strand"),
                                                ("feature_id", None),
                                                ("is_primary_feature", False),
                                            ]
                                        )
                                    ],
                                ),
                                ("feature_collection_name", "unjoined_minus_strand"),
                                ("feature_collection_id", None),
                                ("locus_tag", None),
                                ("feature_collection_type", None),
                                ("sequence_name", "CM021111.1"),
                                ("sequence_guid", None),
                                ("feature_collection_guid", "6332c330-1268-5088-fd99-fc18bc37a9a2"),
                                ("qualifiers", {"name": ["unjoined_minus_strand"]}),
                            ]
                        ),
                    ],
                ),
                (
                    "genes",
                    [
                        OrderedDict(
                            [
                                (
                                    "transcripts",
                                    [
                                        OrderedDict(
                                            [
                                                ("exon_starts", [16174]),
                                                ("exon_ends", [18079]),
                                                ("strand", "MINUS"),
                                                ("cds_starts", None),
                                                ("cds_ends", None),
                                                ("cds_frames", None),
                                                (
                                                    "qualifiers",
                                                    {
                                                        "ncRNA_class": ["other"],
                                                        "locus_tag": ["GI526_G0000001"],
                                                        "product": ["CAT novel prediction: IsoSeq"],
                                                        "note": ["negative strand non-coding gene"],
                                                    },
                                                ),
                                                ("is_primary_tx", False),
                                                ("transcript_id", None),
                                                ("protein_id", None),
                                                ("product", "CAT novel prediction: IsoSeq"),
                                                ("transcript_symbol", None),
                                                ("transcript_type", "ncRNA"),
                                                ("sequence_name", "CM021111.1"),
                                                ("sequence_guid", None),
                                                ("transcript_interval_guid", "7cf2db76-c9db-96e8-c85f-74daf603d9cf"),
                                                ("transcript_guid", None),
                                            ]
                                        )
                                    ],
                                ),
                                ("gene_id", None),
                                ("gene_symbol", None),
                                ("gene_type", "ncRNA"),
                                ("locus_tag", "GI526_G0000001"),
                                ("qualifiers", {"locus_tag": ["GI526_G0000001"]}),
                                ("sequence_name", "CM021111.1"),
                                ("sequence_guid", None),
                                ("gene_guid", "6cb70a2b-c2aa-542d-f382-0a10980cb00a"),
                            ]
                        )
                    ],
                ),
                ("name", "CM021111.1"),
                ("id", None),
                ("sequence_name", "CM021111.1"),
                ("sequence_guid", None),
                ("sequence_path", None),
                (
                    "qualifiers",
                    {
                        "organism": ["Saccharomyces cerevisiae"],
                        "mol_type": ["genomic DNA"],
                        "strain": ["INSC1006"],
                        "db_xref": ["taxon:4932"],
                        "chromosome": ["I"],
                        "country": ["USA: Boulder, CO"],
                    },
                ),
                ("start", 0),
                ("end", 50040),
                ("completely_within", None),
                ("parent_or_seq_chunk_parent", None),
            ]
        )

    def test_caret_coordinates(self, test_data_dir):
        """Handle caret-containing coordinates just fine"""
        genbank = test_data_dir / "caret_coordinates.gbk"
        recs = list(parse_genbank(test_data_dir / genbank))
        c = recs[0].annotation
        assert len(c.genes) == 0
        assert len(c.feature_collections) == 2


class TestSortedParser:
    """Test edge cases for sorted parser"""

    @pytest.mark.parametrize(
        "genbank",
        [
            "INSC1003_missing_gene.gbk",
            "INSC1003_no_genes.gbk",
            "INSC1003_no_genes_misordered.gbk",
        ],
    )
    def test_missing_gene(self, test_data_dir, genbank):
        recs = list(parse_genbank(test_data_dir / genbank, gbk_type=GenBankParserType.SORTED))
        c = recs[0].annotation
        assert len(c.genes) == 8
        assert len([x for x in c.genes if x.gene_type == Biotype.protein_coding]) == 6

    def test_mrna_before_gene(self, test_data_dir):
        genbank = test_data_dir / "INSC1006_chrI_mrna_before_gene.gb"
        recs = list(parse_genbank(test_data_dir / genbank, gbk_type=GenBankParserType.SORTED))
        c = recs[0].annotation
        assert len(c.genes) == 4

    def test_wrong_biotype(self, test_data_dir):
        genbank = test_data_dir / "INSC1006_wrong_feature_type.gbk"
        recs = list(parse_genbank(test_data_dir / genbank, gbk_type=GenBankParserType.SORTED))
        c = recs[0].annotation
        assert len(c.genes) == 1
        gene = c.genes[0]
        assert gene.gene_type == Biotype.ncRNA
        tx = gene.transcripts[0].to_transcript_interval()
        assert tx.is_coding

    def test_misordered(self, test_data_dir):
        genbank = test_data_dir / "INSC1003_misordered.gbk"
        recs = list(parse_genbank(test_data_dir / genbank, gbk_type=GenBankParserType.SORTED))
        c = recs[0].annotation
        assert AnnotationCollectionModel.Schema().dump(c) == OrderedDict(
            [
                ("feature_collections", None),
                (
                    "genes",
                    [
                        OrderedDict(
                            [
                                (
                                    "transcripts",
                                    [
                                        OrderedDict(
                                            [
                                                ("exon_starts", [329]),
                                                ("exon_ends", [2800]),
                                                ("strand", "PLUS"),
                                                ("cds_starts", [334]),
                                                ("cds_ends", [2797]),
                                                ("cds_frames", ["ZERO"]),
                                                ("qualifiers", {"gene": ["thrA"]}),
                                                ("is_primary_tx", False),
                                                ("transcript_id", None),
                                                ("protein_id", None),
                                                ("product", None),
                                                ("transcript_symbol", "thrA"),
                                                ("transcript_type", "protein_coding"),
                                                ("sequence_name", "FEPOIHMA_1"),
                                                ("sequence_guid", None),
                                                ("transcript_interval_guid", "4548a57b-41ec-dcc5-85ca-565412a1d59e"),
                                                ("transcript_guid", None),
                                            ]
                                        )
                                    ],
                                ),
                                ("gene_id", None),
                                ("gene_symbol", "thrA"),
                                ("gene_type", "protein_coding"),
                                ("locus_tag", None),
                                ("qualifiers", {"gene": ["thrA"]}),
                                ("sequence_name", "FEPOIHMA_1"),
                                ("sequence_guid", None),
                                ("gene_guid", "bac88ae4-8de0-147c-a151-887bb370cd3c"),
                            ]
                        ),
                        OrderedDict(
                            [
                                (
                                    "transcripts",
                                    [
                                        OrderedDict(
                                            [
                                                ("exon_starts", [2798]),
                                                ("exon_ends", [3731]),
                                                ("strand", "PLUS"),
                                                ("cds_starts", [2798]),
                                                ("cds_ends", [3731]),
                                                ("cds_frames", ["ZERO"]),
                                                ("qualifiers", {"gene": ["thrB"]}),
                                                ("is_primary_tx", False),
                                                ("transcript_id", None),
                                                ("protein_id", None),
                                                ("product", None),
                                                ("transcript_symbol", "thrB"),
                                                ("transcript_type", "protein_coding"),
                                                ("sequence_name", "FEPOIHMA_1"),
                                                ("sequence_guid", None),
                                                ("transcript_interval_guid", "8a837d19-5fdc-78e4-dd98-8b504d826d46"),
                                                ("transcript_guid", None),
                                            ]
                                        )
                                    ],
                                ),
                                ("gene_id", None),
                                ("gene_symbol", "thrB"),
                                ("gene_type", "protein_coding"),
                                ("locus_tag", None),
                                ("qualifiers", {"gene": ["thrB"]}),
                                ("sequence_name", "FEPOIHMA_1"),
                                ("sequence_guid", None),
                                ("gene_guid", "b3254e70-1495-7b59-d944-b962dc7415eb"),
                            ]
                        ),
                        OrderedDict(
                            [
                                (
                                    "transcripts",
                                    [
                                        OrderedDict(
                                            [
                                                ("exon_starts", [6755]),
                                                ("exon_ends", [6843]),
                                                ("strand", "PLUS"),
                                                ("cds_starts", None),
                                                ("cds_ends", None),
                                                ("cds_frames", None),
                                                ("qualifiers", {"product": ["tRNA-Pro"]}),
                                                ("is_primary_tx", False),
                                                ("transcript_id", None),
                                                ("protein_id", None),
                                                ("product", "tRNA-Pro"),
                                                ("transcript_symbol", None),
                                                ("transcript_type", "tRNA"),
                                                ("sequence_name", "FEPOIHMA_1"),
                                                ("sequence_guid", None),
                                                ("transcript_interval_guid", "54e26268-7e89-5cb1-d9b0-4e9d5944dae6"),
                                                ("transcript_guid", None),
                                            ]
                                        )
                                    ],
                                ),
                                ("gene_id", None),
                                ("gene_symbol", None),
                                ("gene_type", "tRNA"),
                                ("locus_tag", "FEPOIHMA_02080"),
                                ("qualifiers", {"locus_tag": ["FEPOIHMA_02080"]}),
                                ("sequence_name", "FEPOIHMA_1"),
                                ("sequence_guid", None),
                                ("gene_guid", "841a7bc3-79bc-d9ca-9923-5d3f239a2dde"),
                            ]
                        ),
                        OrderedDict(
                            [
                                (
                                    "transcripts",
                                    [
                                        OrderedDict(
                                            [
                                                ("exon_starts", [6899]),
                                                ("exon_ends", [6950]),
                                                ("strand", "PLUS"),
                                                ("cds_starts", None),
                                                ("cds_ends", None),
                                                ("cds_frames", None),
                                                ("qualifiers", {"gene": ["fake"]}),
                                                ("is_primary_tx", False),
                                                ("transcript_id", None),
                                                ("protein_id", None),
                                                ("product", None),
                                                ("transcript_symbol", "fake"),
                                                ("transcript_type", "protein_coding"),
                                                ("sequence_name", "FEPOIHMA_1"),
                                                ("sequence_guid", None),
                                                ("transcript_interval_guid", "f8fba432-72c3-1d17-d327-2f57e5956435"),
                                                ("transcript_guid", None),
                                            ]
                                        )
                                    ],
                                ),
                                ("gene_id", None),
                                ("gene_symbol", "fake"),
                                ("gene_type", "protein_coding"),
                                ("locus_tag", None),
                                ("qualifiers", {"gene": ["fake"]}),
                                ("sequence_name", "FEPOIHMA_1"),
                                ("sequence_guid", None),
                                ("gene_guid", "014d29f5-58e8-b5c1-2b12-d66e7dabf849"),
                            ]
                        ),
                    ],
                ),
                ("name", None),
                ("id", None),
                ("sequence_name", "FEPOIHMA_1"),
                ("sequence_guid", None),
                ("sequence_path", None),
                ("qualifiers", None),
                ("start", 0),
                ("end", 7200),
                ("completely_within", None),
                ("parent_or_seq_chunk_parent", None),
            ]
        )

    def test_misordered_yeast(self, test_data_dir):
        genbank = test_data_dir / "INSC1006_chrI_misordered.gb"
        recs = list(parse_genbank(test_data_dir / genbank, gbk_type=GenBankParserType.SORTED))
        c = recs[0].annotation
        assert AnnotationCollectionModel.Schema().dump(c) == OrderedDict(
            [
                ("feature_collections", None),
                (
                    "genes",
                    [
                        OrderedDict(
                            [
                                (
                                    "transcripts",
                                    [
                                        OrderedDict(
                                            [
                                                ("exon_starts", [16174]),
                                                ("exon_ends", [18079]),
                                                ("strand", "MINUS"),
                                                ("cds_starts", None),
                                                ("cds_ends", None),
                                                ("cds_frames", None),
                                                ("qualifiers", None),
                                                ("is_primary_tx", False),
                                                ("transcript_id", None),
                                                ("protein_id", None),
                                                ("product", None),
                                                ("transcript_symbol", None),
                                                ("transcript_type", "ncRNA"),
                                                ("sequence_name", "CM021111.1"),
                                                ("sequence_guid", None),
                                                ("transcript_interval_guid", "071c9bc1-6eaf-0c54-d52a-586b1fa2c6a5"),
                                                ("transcript_guid", None),
                                            ]
                                        )
                                    ],
                                ),
                                ("gene_id", None),
                                ("gene_symbol", None),
                                ("gene_type", "ncRNA"),
                                ("locus_tag", None),
                                ("qualifiers", {"note": ["gene_after_ncrnA"]}),
                                ("sequence_name", "CM021111.1"),
                                ("sequence_guid", None),
                                ("gene_guid", "e381e2b6-f541-69d2-f00d-0d5339d402a3"),
                            ]
                        ),
                        OrderedDict(
                            [
                                (
                                    "transcripts",
                                    [
                                        OrderedDict(
                                            [
                                                ("exon_starts", [37461]),
                                                ("exon_ends", [39103]),
                                                ("strand", "PLUS"),
                                                ("cds_starts", [37637]),
                                                ("cds_ends", [39011]),
                                                ("cds_frames", ["ZERO"]),
                                                (
                                                    "qualifiers",
                                                    {
                                                        "gene": ["GDH3"],
                                                        "note": ["cds_in_middle", "mrna_before_cds_gene"],
                                                    },
                                                ),
                                                ("is_primary_tx", False),
                                                ("transcript_id", None),
                                                ("protein_id", None),
                                                ("product", None),
                                                ("transcript_symbol", "GDH3"),
                                                ("transcript_type", "protein_coding"),
                                                ("sequence_name", "CM021111.1"),
                                                ("sequence_guid", None),
                                                ("transcript_interval_guid", "76813021-961f-4eaa-e6d8-2ef87b0d90b5"),
                                                ("transcript_guid", None),
                                            ]
                                        )
                                    ],
                                ),
                                ("gene_id", None),
                                ("gene_symbol", "GDH3"),
                                ("gene_type", "protein_coding"),
                                ("locus_tag", None),
                                ("qualifiers", {"gene": ["GDH3"], "note": ["gene_at_end"]}),
                                ("sequence_name", "CM021111.1"),
                                ("sequence_guid", None),
                                ("gene_guid", "151367c7-f12e-172f-f3f8-2ec9c7c0d1f1"),
                            ]
                        ),
                        OrderedDict(
                            [
                                (
                                    "transcripts",
                                    [
                                        OrderedDict(
                                            [
                                                ("exon_starts", [39518]),
                                                ("exon_ends", [40772]),
                                                ("strand", "PLUS"),
                                                ("cds_starts", [39518]),
                                                ("cds_ends", [40772]),
                                                ("cds_frames", ["ZERO"]),
                                                ("qualifiers", {"gene": ["BDH2"], "note": ["cds_in_middle"]}),
                                                ("is_primary_tx", False),
                                                ("transcript_id", None),
                                                ("protein_id", None),
                                                ("product", None),
                                                ("transcript_symbol", "BDH2"),
                                                ("transcript_type", "protein_coding"),
                                                ("sequence_name", "CM021111.1"),
                                                ("sequence_guid", None),
                                                ("transcript_interval_guid", "b8331bbf-c851-ff60-eee4-8e2b6cbbad95"),
                                                ("transcript_guid", None),
                                            ]
                                        )
                                    ],
                                ),
                                ("gene_id", None),
                                ("gene_symbol", "BDH2"),
                                ("gene_type", "protein_coding"),
                                ("locus_tag", None),
                                ("qualifiers", {"gene": ["BDH2"]}),
                                ("sequence_name", "CM021111.1"),
                                ("sequence_guid", None),
                                ("gene_guid", "0529c9ed-3ef4-9642-044c-5498a8e5f2d3"),
                            ]
                        ),
                    ],
                ),
                ("name", None),
                ("id", None),
                ("sequence_name", "CM021111.1"),
                ("sequence_guid", None),
                ("sequence_path", None),
                ("qualifiers", None),
                ("start", 0),
                ("end", 50040),
                ("completely_within", None),
                ("parent_or_seq_chunk_parent", None),
            ]
        )

    def test_overlapping_cds_parser(self, test_data_dir):
        """Overlapping features will be parsed, but the Frames field will be wrong because GenBank files
        do not have the necessary information to properly encode overlapping ORFs."""
        genbank = test_data_dir / "overlapping_cds_feature.gb"
        recs = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(test_data_dir / genbank)))
        assert recs[0].genes[0].transcripts[0].cds.num_blocks > 1
        assert (
            str(recs[0].genes[0].transcripts[0].cds.translate())
            == "MLFAYSGCLAPQCIPDISSFKALPFRDTESRFTTDSSVISSRFSSSFTSSSSKIIIITSIFSSKMDNEHVGASLIVSLSMASLILTNVFSFSSTSYSSQPSDYIACSPSGIDDQPVAEPSGYTPVGSPLHILVVLLLVWMQ*VTNLHLMKANRRMHLPDSSQKFTVILP*LI*FYIYSVSFRSFTYK*APFL"  # noqa: E501
        )


class TestHybridParser:
    def test_hybrid_parser_insc1003(self, test_data_dir):
        genbank = test_data_dir / "INSC1003_hybrid.gbk"
        lt_recs = list(parse_genbank(test_data_dir / genbank, gbk_type=GenBankParserType.LOCUS_TAG))
        hybrid_recs = list(parse_genbank(test_data_dir / genbank, gbk_type=GenBankParserType.HYBRID))
        assert len(hybrid_recs[0].annotation.genes) == 7
        assert len(lt_recs[0].annotation.genes) == 5
        assert len(lt_recs[0].annotation.feature_collections) == 4
        assert hybrid_recs[0].annotation.feature_collections is None

    def test_toy_overlapping_genes_cds_only(self, test_data_dir):
        genbank = test_data_dir / "ToyChr_R64v5_overlapping_genes.gb"
        recs = list(parse_genbank(test_data_dir / genbank, gbk_type=GenBankParserType.HYBRID))
        c = recs[0].annotation
        with open(test_data_dir / "ToyChr_R64v5_overlapping_genes.json") as fh:
            assert AnnotationCollectionModel.Schema().load(json.load(fh)) == c


class TestExceptionsWarnings:
    def test_ambiguous_strand(self, test_data_dir):
        genbank = test_data_dir / "INSC1003_ambiguous_strand.gbk"
        with pytest.warns(StrandViolationWarning):
            recs = list(parse_genbank(test_data_dir / genbank, gbk_type=GenBankParserType.SORTED))
            c = recs[0].annotation
            assert len(c.genes) == 6
        with pytest.warns(StrandViolationWarning):
            recs = list(parse_genbank(test_data_dir / genbank))
            c = recs[0].annotation
            assert len(c.genes) == 6

    @pytest.mark.parametrize(
        "gbk",
        [
            # broken feature
            "broken_coordinates_1.gbk",
            # broken gene
            "broken_coordinates_2.gbk",
        ],
    )
    def test_broken_coordinates(self, test_data_dir, gbk):
        gbk = test_data_dir / gbk
        with pytest.raises(GenBankLocationException):
            with open(gbk, "r") as fh:
                _ = list(parse_genbank(fh))

    def test_duplicate_sequence(self, test_data_dir):
        gbk = test_data_dir / "INSC1006_chrI_duplicate.gbff"
        with pytest.raises(DuplicateSequenceException):
            _ = list(parse_genbank(test_data_dir / gbk))
