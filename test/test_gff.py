from io import StringIO

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from GenBankToLib.gff import GFFRecord, parse_gff, parse_gff3, write_gff, write_gff3


def test_parse_gff3_attributes():
    handle = StringIO(
        "seq1\tfeature\tCDS\t1\t9\t.\t+\t0\t"
        "ID=cds1;product=hypothetical protein;Dbxref=GeneID:1,UniProtKB/Swiss-Prot:P1\n"
        "##FASTA\n"
        ">seq1\n"
        "ATGAAATAA\n"
    )

    records = list(parse_gff3(handle))

    assert records == [
        GFFRecord(
            seqid="seq1",
            source="feature",
            type="CDS",
            start=1,
            end=9,
            score=".",
            strand="+",
            phase="0",
            attributes={
                "ID": ["cds1"],
                "product": ["hypothetical protein"],
                "Dbxref": ["GeneID:1", "UniProtKB/Swiss-Prot:P1"],
            },
        )
    ]


def test_parse_gff2_attributes():
    handle = StringIO("seq1\tfeature\tgene\t1\t9\t.\t+\t.\tgene geneA;locus_tag LT1\n")

    records = list(parse_gff(handle))

    assert records[0].attributes == {
        "gene": ["geneA"],
        "locus_tag": ["LT1"],
    }


def test_write_gff3_records():
    handle = StringIO()
    write_gff(
        [
            GFFRecord(
                seqid="seq1",
                source="feature",
                type="CDS",
                start=1,
                end=9,
                strand="+",
                phase="0",
                attributes={
                    "product": ["hypothetical protein"],
                    "Dbxref": ["GeneID:1", "UniProtKB/Swiss-Prot:P1"],
                },
            )
        ],
        handle,
        gff3=True,
    )

    assert handle.getvalue() == (
        "##gff-version 3\n"
        "seq1\tfeature\tCDS\t1\t9\t.\t+\t0\t"
        "Dbxref=GeneID:1,UniProtKB/Swiss-Prot:P1;product=hypothetical protein\n"
    )


def test_write_gff3_splits_compound_cds_locations():
    seq_record = SeqRecord(Seq("N" * 300), id="seq1")
    seq_record.features = [
        SeqFeature(
            CompoundLocation(
                [
                    SimpleLocation(0, 10, strand=1),
                    SimpleLocation(100, 200, strand=1),
                    SimpleLocation(250, 260, strand=1),
                ]
            ),
            type="CDS",
            qualifiers={"ID": ["cds1"], "codon_start": ["1"]},
        )
    ]
    handle = StringIO()

    write_gff3([seq_record], handle, include_fasta=False)

    records = list(parse_gff3(StringIO(handle.getvalue())))
    assert [(record.start, record.end, record.phase) for record in records] == [
        (1, 10, "0"),
        (101, 200, "2"),
        (251, 260, "1"),
    ]
    assert all(record.type == "CDS" for record in records)
    assert all(record.attributes["ID"] == ["cds1"] for record in records)


def test_write_gff3_splits_origin_spanning_compound_locations():
    seq_record = SeqRecord(Seq("N" * 1000), id="seq1")
    seq_record.features = [
        SeqFeature(
            CompoundLocation(
                [
                    SimpleLocation(899, 1000, strand=1),
                    SimpleLocation(0, 100, strand=1),
                ]
            ),
            type="CDS",
            qualifiers={"ID": ["cds_wrap"], "codon_start": ["2"]},
        )
    ]
    handle = StringIO()

    write_gff3([seq_record], handle, include_fasta=False)

    records = list(parse_gff3(StringIO(handle.getvalue())))
    assert [(record.start, record.end, record.phase) for record in records] == [
        (900, 1000, "1"),
        (1, 100, "2"),
    ]


def test_write_gff3_uses_protein_id_for_compound_feature_id():
    seq_record = SeqRecord(Seq("N" * 30), id="seq1")
    seq_record.features = [
        SeqFeature(
            CompoundLocation(
                [SimpleLocation(0, 9, strand=1), SimpleLocation(20, 30, strand=1)]
            ),
            type="CDS",
            qualifiers={"protein_id": ["protein1"]},
        )
    ]
    handle = StringIO()

    write_gff3([seq_record], handle, include_fasta=False)

    records = list(parse_gff3(StringIO(handle.getvalue())))
    assert [record.attributes["ID"] for record in records] == [
        ["protein1"],
        ["protein1"],
    ]
    assert "ID" not in seq_record.features[0].qualifiers


def test_write_gff3_synthesizes_compound_feature_id():
    seq_record = SeqRecord(Seq("N" * 30), id="seq1")
    seq_record.features = [
        SeqFeature(SimpleLocation(10, 15), type="gene"),
        SeqFeature(
            CompoundLocation([SimpleLocation(0, 5), SimpleLocation(25, 30)]),
            type="misc_feature",
        ),
    ]
    handle = StringIO()

    write_gff3([seq_record], handle, include_fasta=False)

    records = list(parse_gff3(StringIO(handle.getvalue())))
    assert "ID" not in records[0].attributes
    assert [record.attributes["ID"] for record in records[1:]] == [
        ["seq1:misc_feature:2"],
        ["seq1:misc_feature:2"],
    ]
