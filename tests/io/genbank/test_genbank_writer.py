"""
Test GenBank writing.
"""

import pytest

from inscripta.biocantor.io.genbank.parser import parse_genbank, GenBankParserType
from inscripta.biocantor.io.genbank.writer import collection_to_genbank, GenbankFlavor
from inscripta.biocantor.io.parser import ParsedAnnotationRecord
from inscripta.biocantor.gene import AnnotationCollection, GeneInterval, TranscriptInterval
from inscripta.biocantor.location import Strand
from inscripta.biocantor.sequence import Sequence, Parent, Alphabet


@pytest.mark.parametrize(
    "gbk,flavor", [("INSC1006_chrI.gbff", GenbankFlavor.EUKARYOTIC), ("INSC1003.gbk", GenbankFlavor.PROKARYOTIC)]
)
def test_sequences(test_data_dir, tmp_path, gbk, flavor):
    """Parse a genbank, write it to disk, then parse it again and compare."""
    gbk = test_data_dir / gbk
    with open(gbk, "r") as fh:
        collections = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(fh)))

    tmp_gbk = tmp_path / "tmp.gbk"
    with open(tmp_gbk, "w") as fh:
        collection_to_genbank(collections, fh, flavor)

    with open(tmp_gbk, "r") as fh:
        new_collection = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(fh)))

    assert len(collections[0].genes) == len(new_collection[0].genes)

    for gene_a, gene_b in zip(collections[0], new_collection[0]):
        assert gene_a._location == gene_b._location
        tx_a = gene_a.transcripts[0]
        tx_b = gene_b.transcripts[0]
        assert tx_a._location == tx_b._location
        if tx_a.is_coding:
            assert tx_a.cds._location == tx_b.cds._location
            assert tx_a.get_protein_sequence() == tx_b.get_protein_sequence()
        assert tx_a.get_transcript_sequence() == tx_b.get_transcript_sequence()


def test_missing_translation(test_data_dir, tmp_path):
    gbk = test_data_dir / "INSC1003_wrong_missing_translation.gbk"
    with open(gbk, "r") as fh:
        collections = list(
            ParsedAnnotationRecord.parsed_annotation_records_to_model(
                parse_genbank(fh, gbk_type=GenBankParserType.SORTED)
            )
        )

    tmp_gbk = tmp_path / "tmp.gbk"
    with open(tmp_gbk, "w") as fh:
        collection_to_genbank(collections, fh, update_translations=False)

    with open(tmp_gbk, "r") as fh:
        annot = list(
            ParsedAnnotationRecord.parsed_annotation_records_to_model(
                parse_genbank(fh, gbk_type=GenBankParserType.SORTED)
            )
        )[0]
        genes = annot.genes
        assert "translation" not in genes[0].transcripts[0].qualifiers
        assert "translation" in genes[1].transcripts[0].qualifiers
        assert genes[1].transcripts[0].qualifiers["translation"] == {
            "MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGRFADKLPSEPRENIVYQCWERFCQELGK"
            "QIPVAMTLEKNMPIGSGLGSSACSVVAALMAMNEHCGKPLNDTRLLALMGELEGRISGSIHYDNVAPCFLGGMQLMIEE"
            "NDIPELAAKLMKDVIAEPYRERLLPGFRQARQAVAEIGAVASGISGSGPTLFALCDKPDTAQRVADWLGKNYLQNQEGF"
            "VHICRLDTAGARVLEN"
        }

    # now export to file and force the translations to be recalculated
    tmp_gbk = tmp_path / "tmp.gbk"
    with open(tmp_gbk, "w") as fh:
        collection_to_genbank(collections, fh, update_translations=True)

    with open(tmp_gbk, "r") as fh:
        annot = list(
            ParsedAnnotationRecord.parsed_annotation_records_to_model(
                parse_genbank(fh, gbk_type=GenBankParserType.SORTED)
            )
        )[0]
        genes = annot.genes
        assert genes[0].transcripts[0].qualifiers["translation"] == {
            "MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDALPNISDAERIFAELLTGLAAA"
            "QPGFPLAQLKTFVDQEFAQIKHVLHGISLLGQCPDSINAALICRGEKMSIAIMAGVLEARGHNVTVIDPVEKLLAVGHYLE"
            "STVDIAESTRRIAASRIPADHMVLMAGFTAGNEKGELVVLGRNGSDYSAAVLAACLRADCCEIWTDVDGVYTCDPRQVPDAR"
            "LLKSMSYQEAMELSYFGAKVLHPRTITPIAQFQIPCLIKNTGNPQAPGTLIGASRDEDELPVKGISNLNNMAMFSVSGPGMK"
            "GMVGMAARVFAAMSRARISVVLITQSSSEYSISFCVPQSDCVRAERAMQEEFYLELKEGLLEPLAVTERLAIISVVGDGMRT"
            "LRGISAKFFAALARANINIVAIAQGSSERSISVVVNNDDATTGVRVTHQMLFNTDQVIEVFVIGVGGVGGALLEQLKRQQSWL"
            "KNKHIDLRVCGVANSKALLTNVHGLNLENWQEELAQAKEPFNLGRLIRLVKEYHLLNPVIVDCTSSQAVADQYADFLREGFHV"
            "VTPNKKANTSSMDYYHLLRHAAEKSRRKFLYDTNVGAGLPVIENLQNLLNAGDELMKFSGILSGSLSYIFGKLDEGMSFSEATT"
            "LAREMGYTEPDPRDDLSGMDVARKLLILARETGRELELADIEIEPVLPAEFNAEGDVAAFMANLSQLDDLFAARVAKARDEGKVL"
            "RYVGNIDEDGACRVKIAEVDGNDPLFKVKNGENALAFYSHYYQPLPLVLRGYGAGNDVTAAGVFADLLRTLSWKLGV*"
        }
        assert genes[1].transcripts[0].qualifiers["translation"] == {
            "MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGRFADKLPSEPRENIVYQCWERFCQELGK"
            "QIPVAMTLEKNMPIGSGLGSSACSVVAALMAMNEHCGKPLNDTRLLALMGELEGRISGSIHYDNVAPCFLGGMQLMIEE"
            "NDIISQQVPGFDEWLWVLAYPGIKVSTAEARAILPAQYRRQDCIAHGRHLAGFIHACYSRQPELAAKLMKDVIAEPYRE"
            "RLLPGFRQARQAVAEIGAVASGISGSGPTLFALCDKPDTAQRVADWLGKNYLQNQEGFVHICRLDTAGARVLEN*"
        }


def test_locus_tag(tmp_path):
    """Null locus_tag should lead to gene symbol as locus tag on every feature"""
    parent = Parent(sequence=Sequence("A" * 10, Alphabet.NT_EXTENDED_GAPPED, id="chr"))
    a = AnnotationCollection(
        None,
        [
            GeneInterval(
                [TranscriptInterval([0], [10], Strand.PLUS, parent_or_seq_chunk_parent=parent)],
                parent_or_seq_chunk_parent=parent,
                gene_symbol="mygene",
            )
        ],
        parent_or_seq_chunk_parent=parent,
        sequence_name="chr",
    )
    tmp_gb = tmp_path / "tmp.gb"
    collection_to_genbank([a], tmp_gb, genbank_type=GenbankFlavor.EUKARYOTIC)
    assert (
        open(tmp_gb).read()
        == """LOCUS       chr                       10 bp    DNA              UNK 01-JAN-1980
DEFINITION  GenBank produced by BioCantor.
ACCESSION   chr
VERSION     chr
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     gene            1..10
                     /gene_name="mygene"
                     /gene_biotype="unspecified"
                     /gene="mygene"
                     /locus_tag="mygene"
     misc_RNA        1..10
                     /transcript_biotype="unspecified"
                     /gene="mygene"
                     /locus_tag="mygene"
ORIGIN
        1 aaaaaaaaaa
//
"""
    )
