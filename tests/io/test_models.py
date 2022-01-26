"""
Test importing and exporting AnnotationCollectionModel objects from JSON, with and without sequence info.

All of these tests are based on the INSC1003.gbk file.
"""
import pytest
import json
from uuid import UUID
from inscripta.biocantor.io.models import (
    AnnotationCollectionModel,
    TranscriptIntervalModel,
    GeneIntervalModel,
    Strand,
    CDSFrame,
    Biotype,
    SequenceType,
    Alphabet,
    ParentModel,
)


model_without_parent = AnnotationCollectionModel(
    feature_collections=[],
    genes=[
        GeneIntervalModel(
            transcripts=[
                TranscriptIntervalModel(
                    exon_starts=[334],
                    exon_ends=[2797],
                    strand=Strand.PLUS,
                    cds_starts=[334],
                    cds_ends=[2797],
                    cds_frames=[CDSFrame.ZERO],
                    qualifiers=None,
                    is_primary_tx=False,
                    transcript_id=None,
                    protein_id="X:FEPOIHMA_00001",
                    product="aspartate kinase / homoserine dehydrogenase",
                    transcript_symbol="thrA",
                    transcript_type=Biotype.protein_coding,
                    sequence_name="FEPOIHMA_1",
                    sequence_guid=None,
                    transcript_interval_guid=UUID("6fd4a86d-62b0-eafd-70f5-33cb2d575769"),
                    transcript_guid=None,
                )
            ],
            gene_id=None,
            gene_symbol="thrA",
            gene_type=Biotype.protein_coding,
            locus_tag="FEPOIHMA_00001",
            qualifiers=None,
            sequence_name="FEPOIHMA_1",
            sequence_guid=None,
            gene_guid=UUID("ff1fa9c6-fc20-f699-7164-f0c2e235a487"),
        ),
        GeneIntervalModel(
            transcripts=[
                TranscriptIntervalModel(
                    exon_starts=[2798],
                    exon_ends=[3731],
                    strand=Strand.PLUS,
                    cds_starts=[2798],
                    cds_ends=[3731],
                    cds_frames=[CDSFrame.ZERO],
                    qualifiers=None,
                    is_primary_tx=False,
                    transcript_id=None,
                    protein_id="X:FEPOIHMA_00002",
                    product="homoserine kinase",
                    transcript_symbol="thrB",
                    transcript_type=Biotype.protein_coding,
                    sequence_name="FEPOIHMA_1",
                    sequence_guid=None,
                    transcript_interval_guid=UUID("a51bd755-1683-7c7a-c9fb-1ae7f49cb249"),
                    transcript_guid=None,
                )
            ],
            gene_id=None,
            gene_symbol="thrB",
            gene_type=Biotype.protein_coding,
            locus_tag="FEPOIHMA_00002",
            qualifiers={"gene": ["thrB"], "locus_tag": ["FEPOIHMA_00002"]},
            sequence_name="FEPOIHMA_1",
            sequence_guid=None,
            gene_guid=UUID("f497ba2e-0b43-dae5-5300-56a7976de0c8"),
        ),
        GeneIntervalModel(
            transcripts=[
                TranscriptIntervalModel(
                    exon_starts=[3731],
                    exon_ends=[5018],
                    strand=Strand.PLUS,
                    cds_starts=[3731],
                    cds_ends=[5018],
                    cds_frames=[CDSFrame.ZERO],
                    qualifiers=None,
                    is_primary_tx=False,
                    transcript_id=None,
                    protein_id="X:FEPOIHMA_00003",
                    product="threonine synthase",
                    transcript_symbol="thrC",
                    transcript_type=Biotype.protein_coding,
                    sequence_name="FEPOIHMA_1",
                    sequence_guid=None,
                    transcript_interval_guid=UUID("3008ee04-fc1f-9a03-47f1-092e9b527c92"),
                    transcript_guid=None,
                )
            ],
            gene_id=None,
            gene_symbol="thrC",
            gene_type=Biotype.protein_coding,
            locus_tag="FEPOIHMA_00003",
            qualifiers={"gene": ["thrC"], "locus_tag": ["FEPOIHMA_00003"]},
            sequence_name="FEPOIHMA_1",
            sequence_guid=None,
            gene_guid=UUID("329eac8f-31b9-4bc8-20b9-ff1210639449"),
        ),
        GeneIntervalModel(
            transcripts=[
                TranscriptIntervalModel(
                    exon_starts=[5230],
                    exon_ends=[5527],
                    strand=Strand.PLUS,
                    cds_starts=[5230],
                    cds_ends=[5527],
                    cds_frames=[CDSFrame.ZERO],
                    qualifiers=None,
                    is_primary_tx=False,
                    transcript_id=None,
                    protein_id="X:FEPOIHMA_00004",
                    product="hypothetical protein",
                    transcript_symbol="yaaX",
                    transcript_type=Biotype.protein_coding,
                    sequence_name="FEPOIHMA_1",
                    sequence_guid=None,
                    transcript_interval_guid=UUID("078847e2-d770-83e0-e0ff-928bc7d5025f"),
                    transcript_guid=None,
                )
            ],
            gene_id=None,
            gene_symbol="yaaX",
            gene_type=Biotype.protein_coding,
            locus_tag="FEPOIHMA_00004",
            qualifiers={"gene": ["yaaX"], "locus_tag": ["FEPOIHMA_00004"]},
            sequence_name="FEPOIHMA_1",
            sequence_guid=None,
            gene_guid=UUID("8b2e3283-604b-0233-c770-115889d84bbe"),
        ),
        GeneIntervalModel(
            transcripts=[
                TranscriptIntervalModel(
                    exon_starts=[5679],
                    exon_ends=[6456],
                    strand=Strand.MINUS,
                    cds_starts=[5679],
                    cds_ends=[6456],
                    cds_frames=[CDSFrame.ZERO],
                    qualifiers=None,
                    is_primary_tx=False,
                    transcript_id=None,
                    protein_id="X:FEPOIHMA_00005",
                    product="protein that reduces intracellular iron levels under peroxide stress",
                    transcript_symbol="yaaA",
                    transcript_type=Biotype.protein_coding,
                    sequence_name="FEPOIHMA_1",
                    sequence_guid=None,
                    transcript_interval_guid=UUID("6e1ecab2-29fe-7c9c-8394-30e55172d989"),
                    transcript_guid=None,
                )
            ],
            gene_id=None,
            gene_symbol="yaaA",
            gene_type=Biotype.protein_coding,
            locus_tag="FEPOIHMA_00005",
            qualifiers={"gene": ["yaaA"], "locus_tag": ["FEPOIHMA_00005"]},
            sequence_name="FEPOIHMA_1",
            sequence_guid=None,
            gene_guid=UUID("93599885-f0c0-9818-ae68-64816e366a77"),
        ),
        GeneIntervalModel(
            transcripts=[
                TranscriptIntervalModel(
                    exon_starts=[6755],
                    exon_ends=[6843],
                    strand=Strand.PLUS,
                    cds_starts=None,
                    cds_ends=None,
                    cds_frames=None,
                    qualifiers=None,
                    is_primary_tx=False,
                    transcript_id=None,
                    protein_id=None,
                    product="tRNA-Pro",
                    transcript_symbol=None,
                    transcript_type=Biotype.tRNA,
                    sequence_name="FEPOIHMA_1",
                    sequence_guid=None,
                    transcript_interval_guid=UUID("6d44fc35-9f6d-9bb2-c0d9-80fff74ee031"),
                    transcript_guid=None,
                )
            ],
            gene_id=None,
            gene_symbol=None,
            gene_type=Biotype.tRNA,
            locus_tag="FEPOIHMA_02080",
            qualifiers={"locus_tag": ["FEPOIHMA_02080"]},
            sequence_name="FEPOIHMA_1",
            sequence_guid=None,
            gene_guid=UUID("1adc918e-fcb2-18ff-384b-e7d2a109120f"),
        ),
        GeneIntervalModel(
            transcripts=[
                TranscriptIntervalModel(
                    exon_starts=[6955],
                    exon_ends=[7066],
                    strand=Strand.MINUS,
                    cds_starts=None,
                    cds_ends=None,
                    cds_frames=None,
                    qualifiers=None,
                    is_primary_tx=False,
                    transcript_id=None,
                    protein_id=None,
                    product="5S ribosomal RNA",
                    transcript_symbol=None,
                    transcript_type=Biotype.rRNA,
                    sequence_name="FEPOIHMA_1",
                    sequence_guid=None,
                    transcript_interval_guid=UUID("9a3d0d0d-3edd-ba8b-51c9-78377ffa7131"),
                    transcript_guid=None,
                )
            ],
            gene_id=None,
            gene_symbol=None,
            gene_type=Biotype.rRNA,
            locus_tag="FEPOIHMA_02457",
            qualifiers=None,
            sequence_name="FEPOIHMA_1",
            sequence_guid=None,
            gene_guid=UUID("b8a03f84-ba9d-6209-88e4-38f1d54f7e65"),
        ),
    ],
    name="FEPOIHMA_1",
    id=None,
    sequence_name="FEPOIHMA_1",
    sequence_guid=None,
    sequence_path=None,
    qualifiers={"organism": ["Genus species"], "mol_type": ["genomic DNA"], "strain": ["strain"]},
    start=0,
    end=7200,
    completely_within=None,
    parent_or_seq_chunk_parent=None,
)

parent_model = ParentModel(
    seq="GCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGGGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTTGACGGGACTCGCCGCCGCCCAGCCGGGATTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACGCTGCGCTGATTTGCCGTGGCGAGAAAATGTCGATCGCCATTATGGCCGGCGTATTAGAAGCGCGCGGTCACAACGTTACCGTTATCGATCCGGTCGAAAAACTGCTGGCAGTGGGGCATTACCTCGAATCTACCGTCGATATTGCTGAGTCCACCCGCCGTATTGCGGCAAGTCGCATTCCGGCTGATCACATGGTGCTGATGGCAGGTTTCACCGCCGGTAATGAAAAAGGCGAACTGGTGGTACTTGGACGCAACGGTTCCGACTACTCCGCGGCGGTGCTGGCTGCCTGTTTACGCGCCGATTGTTGCGAGATTTGGACGGACGTTGACGGGGTCTATACCTGCGACCCGCGTCAGGTGCCCGATGCGAGGTTGTTGAAGTCGATGTCCTACCAGGAAGCGATGGAGCTTTCCTACTTCGGCGCTAAAGTTCTTCACCCCCGCACCATTACCCCCATCGCCCAGTTCCAGATCCCTTGCCTGATTAAAAATACCGGAAATCCTCAAGCTCCAGGTACGCTCATTGGTGCCAGCCGTGATGAAGACGAATTACCGGTCAAGGGCATTTCCAATCTGAATAATATGGCAATGTTCAGCGTTTCCGGCCCGGGGATGAAAGGGATGGTTGGCATGGCGGCGCGCGTGTTTGCAGCGATGTCACGCGCCCGTATTTCCGTGGTGCTGATTACGCAATCATCTTCCGAATACAGTATCAGTTTCTGCGTTCCGCAAAGCGACTGTGTGCGAGCTGAACGGGCAATGCAGGAAGAGTTCTACCTGGAACTGAAAGAAGGCTTACTGGAGCCGCTGGCGGTGACGGAACGGCTGGCCATTATCTCGGTGGTAGGTGATGGTATGCGCACCTTGCGTGGGATCTCGGCGAAATTCTTTGCCGCGCTGGCCCGCGCCAATATCAACATTGTCGCCATTGCTCAGGGATCTTCTGAACGCTCAATCTCTGTCGTGGTAAATAACGATGATGCGACCACTGGCGTGCGCGTTACTCATCAGATGCTGTTCAATACCGATCAGGTTATCGAAGTGTTTGTGATTGGCGTCGGTGGCGTTGGCGGTGCGCTGCTGGAGCAACTGAAGCGTCAACAAAGCTGGCTGAAGAATAAACATATCGACTTACGTGTCTGCGGTGTTGCCAACTCGAAGGCACTGCTCACCAATGTGCATGGCCTAAATCTGGAAAACTGGCAGGAAGAACTGGCGCAAGCCAAAGAGCCGTTTAATCTCGGGCGCTTAATTCGCCTCGTGAAAGAATATCATCTGCTGAACCCGGTCATTGTTGACTGCACTTCCAGCCAGGCAGTGGCGGATCAATATGCCGACTTCTTGCGCGAAGGTTTCCACGTTGTCACGCCGAACAAAAAGGCCAACACCTCGTCGATGGATTACTACCATCTGTTGCGTCATGCGGCGGAAAAATCGCGGCGTAAATTCCTCTATGACACCAACGTTGGGGCTGGATTACCGGTTATTGAGAACCTGCAAAATCTGCTCAATGCTGGTGATGAATTGATGAAGTTCTCCGGCATTCTTTCAGGTTCGCTTTCTTATATCTTCGGCAAGTTAGACGAAGGCATGAGTTTCTCCGAGGCGACTACTCTGGCGCGGGAAATGGGTTATACCGAACCGGATCCGCGAGATGATCTTTCTGGTATGGATGTAGCGCGTAAGCTATTGATTCTCGCTCGTGAAACGGGACGTGAACTGGAGCTGGCGGATATTGAAATTGAACCTGTGCTGCCCGCAGAGTTTAACGCTGAGGGTGATGTTGCCGCTTTTATGGCGAATCTGTCACAGCTCGACGATCTCTTTGCCGCGCGCGTGGCGAAGGCCCGTGATGAAGGAAAAGTTTTGCGCTATGTTGGCAATATTGATGAAGATGGTGCCTGCCGCGTGAAGATTGCCGAAGTGGATGGTAATGATCCGCTGTTCAAAGTGAAAAATGGCGAAAACGCCCTGGCCTTTTATAGCCACTATTATCAGCCGCTGCCGTTGGTGCTGCGCGGATATGGTGCGGGCAATGACGTTACAGCTGCCGGTGTCTTTGCCGATCTGCTACGTACCCTCTCATGGAAGTTAGGAGTCTGACATGGTTAAAGTTTATGCCCCGGCTTCCAGTGCCAATATGAGCGTCGGGTTTGATGTGCTCGGGGCGGCGGTGACACCTGTTGATGGTGCATTGCTCGGAGATGTAGTCACGGTTGAGGCGGCAGAGACATTCAGTCTCAACAACCTCGGACGCTTTGCCGATAAGCTGCCGTCAGAACCACGGGAAAATATCGTTTATCAGTGCTGGGAGCGTTTTTGCCAGGAGCTTGGCAAGCAAATTCCAGTGGCGATGACTCTGGAAAAGAATATGCCAATCGGTTCGGGCTTAGGCTCCAGCGCCTGTTCGGTGGTCGCGGCGCTGATGGCGATGAATGAACACTGTGGCAAGCCGCTTAATGACACTCGTTTGCTGGCTTTGATGGGCGAGCTGGAAGGACGAATCTCCGGCAGCATTCATTACGACAACGTGGCACCGTGTTTTCTTGGTGGTATGCAGTTGATGATCGAAGAAAACGACATCATCAGCCAGCAAGTGCCAGGGTTTGATGAGTGGCTGTGGGTGCTGGCGTATCCGGGGATTAAAGTCTCGACGGCAGAAGCCAGGGCTATTTTACCGGCGCAGTATCGCCGCCAGGATTGCATTGCGCACGGGCGACATCTGGCTGGCTTCATTCACGCCTGCTATTCCCGTCAGCCTGAGCTTGCCGCGAAGCTGATGAAAGATGTTATCGCTGAACCCTACCGTGAACGGTTACTGCCTGGCTTCCGGCAGGCGCGGCAGGCGGTTGCGGAAATCGGCGCGGTAGCGAGCGGTATCTCCGGCTCCGGCCCGACCTTGTTCGCTCTGTGTGACAAGCCGGATACCGCCCAGCGCGTTGCCGACTGGTTGGGTAAGAACTACCTGCAAAATCAGGAAGGTTTTGTTCATATTTGCCGGCTGGATACGGCGGGCGCACGAGTACTGGAAAACTAAATGAAACTCTACAATCTGAAAGATCACAATGAGCAGGTCAGCTTTGCGCAAGCCGTAACCCAGGGGTTGGGCAAAAATCAGGGGCTGTTTTTTCCGCACGACCTGCCGGAATTCAGCCTGACTGAAATTGATGAGATGCTGAAGCTGGATTTTGTCACCCGCAGTGCGAAGATCCTCTCGGCGTTTATTGGTGATGAAATCCCGCAGGAAATCCTGGAAGAGCGCGTGCGCGCGGCGTTTGCCTTCCCGGCTCCGGTCGCCAATGTTGAAAGCGATGTCGGTTGTCTGGAATTGTTCCACGGGCCAACGCTGGCATTTAAAGATTTCGGCGGTCGCTTTATGGCACAAATGCTGACCCATATTGCGGGCGATAAGCCAGTGACCATTCTGACCGCGACCTCCGGTGATACCGGAGCGGCAGTGGCTCATGCTTTCTACGGTTTACCGAATGTGAAAGTGGTTATCCTCTATCCACGAGGCAAAATCAGTCCACTGCAAGAAAAACTGTTCTGTACATTGGGCGGCAATATCGAAACTGTTGCCATCGACGGCGATTTCGATGCCTGTCAGGCGCTGGTGAAGCAGGCGTTTGATGATGAAGAGCTGAAAGTGGCGCTGGGGTTAAACTCAGCTAACTCGATTAACATCAGCCGTTTGCTGGCGCAGATTTGCTACTACTTTGAAGCAGTTGCGCAGCTGCCGCAGGAAGCGCGCAACCAGCTGGTTGTCTCGGTGCCAAGCGGAAACTTCGGCGATTTGACGGCGGGTCTGCTGGCGAAGTCACTCGGTCTGCCGGTGAAACGTTTTATTGCTGCGACCAACGTGAACGATACCGTGCCACGTTTCCTGCACGACGGTCAGTGGTCACCCAAAGCGACTCAGGCGACGTTATCCAACGCGATGGACGTGAGTCAGCCGAACAACTGGCCGCGTGTGGAAGAGTTGTTCCGCCGCAAAATCTGGCAACTGAAAGAGCTGGGTTATGCAGCCGTGGATGATGAAACCACGCAACAGACAATGCGTGAGTTAAAAGAACTGGGCTACACCTCGGAGCCGCACGCTGCCGTAGCGTATCGTGCGCTGCGTGACCAGTTGAATCCAGGCGAATATGGCTTGTTCCTCGGCACCGCGCATCCGGCGAAATTTAAAGAGAGCGTGGAAGCGATTCTCGGTGAAACGTTGGATCTGCCAAAAGAGCTGGCAGAACGTGCTGATTTACCCTTGCTTTCACATAATCTGCCCGCCGATTTTGCTGCGTTGCGTAAATTGATGATGAATCATCAGTAAAATCTATTCATTATCTCAATCAGGCCGGGTTTGCTTTTATGCAGCCGGCTTTTTTATGAAGAAATTATGGAGAAAAACGACAGGGAAAAAGGAGAAATTCTCAATAAATGCGGTAACTTAGAGATTAGGATTGCGGAGAATAACAACCGTCGTTCTCATCGCGTAATCTCCGGATATCGACCCATAACGGGCAATGATAAAAGGAGTAACCTATGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGATAATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACATTATGAATGGCGAGGCAATCGCTGGCACCCACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAAATGACAAATGCCGGGTAACAATCCGGCATTCAGCGCCTGATGCGACGCTGGCGCGTCTTATCAGGCCTACGTGAATTCTGCAATATATTGAATCTGCATGCTTTTGTAGGCCGGATAAGGCGTTCACGCCGCATCCGGCATTGACTGCAAACTTAACGCTGCTCGTAGCGTTTAAACACCAGTTCGCCATTGCTGGAGGAAGCTTCATCAAAGAAGTAACCTTCGCTATTAAAACCAGTCAGTTGCTCTGGTTTGGTCAGCCGATTTTCAATAATAAAACGACTCATCAGACCGCGTGCTTTCTTAGCGTAGAAGCTGATTATCTTAAATTTGCCGTTCTTCTCATCGAGGAACACCGGCTTGATAATCTCGGCATTCAATTTCTTCGGCTTCACCGATTTAAAATACTCATCTGACGCCAGATTAATCACCACATTATCGCCTTGTGCTGCGAGCGCCTCGTTCAGCTTGTTGGTGATGATATCTCCCCAGAATTGATACAGATCTTTCCCTCGGGCATTCTCAAGACGGATCCCCATTTCCAGACGATAAGGCTGCATTAAATCGAGCGGGCGCAGTACGCCATACAAGCCGGAAAGCATTCGCAAATGCTGTTGGGCAAAATCGAAATCGTCTTCGCTGAAGGTTTCGGCCTGCAAGCCGGTGTAGACATCACCTTTAAACGCCAGAATCGCCTGGCGGGCATTCTCCGGCGTGAAATCTGGCTGCCAGTCATGAAAGCGAGCGGCGTTGATACCCGCCAGTTTGTCGCTGATGCGCATCAGCGTGCTAATCTGCGGAGGCGTCAGTTTCCGCGCTTCATGGATCAACTGCTGGGAATTGTCTAACAGCTCCGGCAGCGTATAGCGCGTGGTGGTCAACGGGCTTTGGTAATCAAGCGTTTTCGCAGGTGAAATAAGAATCAGCATATCCAGTCCTTGCAGGAAATTTATGCCGACTTTAGCAAAAAAAGAGAATGAGTTGATCGATAGTTGTGATTACTCCTGCGAAACATCATCCCACGCGTCCGGAGAAAGCTGGCGGCCGATATCCGGATAACGCAACGGATCAAACACCGGGCGCACGCCGAGTTTACGCTGGCGTAGATAATCACTGGCAATGGTATGAACCACAGGCGAGAGCAGTAAAATGGCGGTCAAATTGGTAATAGCCATGCAGGCCATTATGATATCTGCCAGTTGCCACATCAGCGGAAGACTTAGCAAGGTGCCGCCGATGACCGTTGCGAAGGTGCAGATCCGCAAACACCAGATCGCTTTAGGGTTGTTCAGGCGTAAAAAGAAGAGATTGTTTTCGGCGTAAATGTAGTTGGCAACGATGGAGCTGAAGGCAAACAGAATAACCACGAGGGTAACAAACTCAGCACCCCAGGAACCCATTAACACCCGCATCGCCTTCTGGATAAGCTGAATACCTTCCAGCGGCATGTAGGTTGTGCCGTTACCCGCCAGTAATATCAGCATGGCGCTTGCCGTACAGATGACCAGGGTGTCGATAAAAATGCCAATCATCTGGACAATCCCTTGCGCTGCCGGATGCGGAGGCCAGGACGCCGCTGCCGCTGCCGCGTTTGGCGTCGACCCCATTCCCGCCTCATTGGAAAACATACTGCGCTGAAAACCGTTAGTAATCGCCTGGCTTAAGGTATATCC",  # noqa: E501
    alphabet=Alphabet.NT_EXTENDED_GAPPED,
    sequence_name="FEPOIHMA_1",
    type=SequenceType.CHROMOSOME,
    start=0,
    end=7200,
    strand=Strand.PLUS,
)

dumped = AnnotationCollectionModel.Schema().dump(model_without_parent)
dumped["parent_or_seq_chunk_parent"] = ParentModel.Schema().dump(parent_model)
model_with_parent = AnnotationCollectionModel.Schema().load(dumped)


class TestAnnotationCollectionModel:
    @pytest.mark.parametrize(
        "json_file,model",
        [
            (
                "INSC1003.json",
                model_without_parent,
            ),
            ("INSC1003_with_parent.json", model_with_parent),
        ],
    )
    def test_export(self, test_data_dir, json_file, model):
        """JSON serialized previously is the same value as dumping a model to disk"""
        with open(test_data_dir / json_file) as fh:
            ref = json.load(fh)
        assert ref == AnnotationCollectionModel.Schema().dump(model)

    @pytest.mark.parametrize(
        "json_file,model",
        [
            (
                "INSC1003.json",
                model_without_parent,
            ),
            ("INSC1003_with_parent.json", model_with_parent),
        ],
    )
    def test_import(self, test_data_dir, json_file, model):
        """Converting the model to AnnotationCollection is the same as loading from disk then converting"""
        ac = model.to_annotation_collection()
        with open(test_data_dir / json_file) as fh:
            serialized_ac_model = AnnotationCollectionModel.Schema().loads(fh.read())
        serialized_ac = serialized_ac_model.to_annotation_collection()
        assert serialized_ac == ac
        assert model == serialized_ac_model

    def test_from_annotation_collection_without_parent(self):
        """Convert an AnnotationCollection back to an AnnotationCollectionModel"""
        ac = model_without_parent.to_annotation_collection()
        new_model = AnnotationCollectionModel.from_annotation_collection(ac)
        assert new_model == model_without_parent
        # exporting parent does nothing because it is null
        new_model_with_parent = AnnotationCollectionModel.from_annotation_collection(ac, export_parent=True)
        assert new_model_with_parent == model_without_parent

    def test_from_annotation_collection_with_parent(self):
        """Convert an AnnotationCollection back to an AnnotationCollectionModel"""
        ac = model_with_parent.to_annotation_collection()
        new_model = AnnotationCollectionModel.from_annotation_collection(ac, export_parent=True)
        assert new_model == model_with_parent
        # not exporting parent leads to it looking like the non-parent version
        new_model_without_parent = AnnotationCollectionModel.from_annotation_collection(ac, export_parent=False)
        assert new_model_without_parent == model_without_parent
