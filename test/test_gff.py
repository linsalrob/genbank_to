from io import StringIO

from GenBankToLib.gff import GFFRecord, parse_gff, parse_gff3, write_gff


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
