"""
Small GFF/GFF3 reader and writer utilities.
"""
from dataclasses import dataclass, field
from typing import Iterable, Iterator, TextIO
from urllib.parse import quote, unquote


@dataclass
class GFFRecord:
    """A single GFF feature row."""

    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: str = "."
    strand: str = "."
    phase: str = "."
    attributes: dict[str, list[str]] = field(default_factory=dict)


def parse_attributes(value: str, gff3: bool = True) -> dict[str, list[str]]:
    """
    Parse a GFF or GFF3 attribute column.

    GFF3 uses key=value pairs. Older GFF files often use key value pairs, so this
    accepts both forms.
    """

    if value == ".":
        return {}

    attributes: dict[str, list[str]] = {}
    for item in value.rstrip().split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, raw_value = item.split("=", 1)
        elif " " in item and not gff3:
            key, raw_value = item.split(" ", 1)
        else:
            key, raw_value = item, ""
        attributes[unquote(key)] = [unquote(v) for v in raw_value.split(",")]
    return attributes


def parse_gff(handle: TextIO, gff3: bool = False) -> Iterator[GFFRecord]:
    """Yield feature rows from a GFF or GFF3 file handle."""

    for line in handle:
        if line.startswith("##FASTA"):
            break
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 9:
            raise ValueError(f"Expected 9 GFF columns, found {len(parts)}: {line.rstrip()}")
        yield GFFRecord(
            seqid=parts[0],
            source=parts[1],
            type=parts[2],
            start=int(parts[3]),
            end=int(parts[4]),
            score=parts[5],
            strand=parts[6],
            phase=parts[7],
            attributes=parse_attributes(parts[8], gff3=gff3),
        )


def parse_gff3(handle: TextIO) -> Iterator[GFFRecord]:
    """Yield feature rows from a GFF3 file handle."""

    yield from parse_gff(handle, gff3=True)


def format_attributes(attributes: dict[str, Iterable[object]], gff3: bool = True) -> str:
    """Format attributes for a GFF or GFF3 row."""

    if not attributes:
        return "."

    items = []
    for key in sorted(attributes):
        raw_values = attributes[key]
        values = list(raw_values)
        if not values:
            values = [""]
        encoded_values = ",".join(_escape(value) for value in values)
        if gff3:
            items.append(f"{_escape(key)}={encoded_values}")
        else:
            items.append(f"{_escape(key)} {encoded_values}")
    return ";".join(items)


def write_gff(records: Iterable[GFFRecord], handle: TextIO, gff3: bool = False) -> None:
    """Write plain GFFRecord objects as GFF or GFF3."""

    if gff3:
        handle.write("##gff-version 3\n")
    for record in records:
        handle.write(_format_record(record, gff3=gff3))


def write_gff3(seq_records: Iterable[object], handle: TextIO, include_fasta: bool = True) -> None:
    """
    Write Biopython SeqRecord objects as GFF3.

    This intentionally covers the GenBank-to-GFF path used by this package:
    record annotations, sequence-region pragmas, feature qualifiers, and an
    optional FASTA trailer.
    """

    records = list(seq_records)
    handle.write("##gff-version 3\n")
    for seq_record in records:
        handle.write(f"##sequence-region {seq_record.id} 1 {len(seq_record.seq)}\n")
        annotation_attributes = _annotation_attributes(seq_record)
        if annotation_attributes:
            annotation_record = GFFRecord(
                seqid=seq_record.id,
                source="annotation",
                type="remark",
                start=1,
                end=len(seq_record.seq),
                attributes=annotation_attributes,
            )
            handle.write(_format_record(annotation_record, gff3=True))
        for feature in seq_record.features:
            handle.write(_format_feature(seq_record, feature))

    if include_fasta and records:
        handle.write("##FASTA\n")
        for seq_record in records:
            description = getattr(seq_record, "description", "") or seq_record.id
            handle.write(f">{seq_record.id} {description}\n")
            sequence = str(seq_record.seq)
            for index in range(0, len(sequence), 60):
                handle.write(f"{sequence[index:index + 60]}\n")


def _format_record(record: GFFRecord, gff3: bool = True) -> str:
    return "\t".join(
        [
            record.seqid,
            record.source,
            record.type,
            str(record.start),
            str(record.end),
            record.score,
            record.strand,
            record.phase,
            format_attributes(record.attributes, gff3=gff3),
        ]
    ) + "\n"


def _format_feature(seq_record: object, feature: object) -> str:
    start = int(feature.location.start) + 1
    end = int(feature.location.end)
    strand = _strand(feature)
    phase = "."
    if feature.type == "CDS":
        codon_start = feature.qualifiers.get("codon_start", ["1"])[0]
        phase = str(int(codon_start) - 1)
    record = GFFRecord(
        seqid=seq_record.id,
        source="feature",
        type=feature.type,
        start=start,
        end=end,
        strand=strand,
        phase=phase,
        attributes=feature.qualifiers,
    )
    return _format_record(record, gff3=True)


def _annotation_attributes(seq_record: object) -> dict[str, list[object]]:
    attributes: dict[str, list[object]] = {}
    for key, value in seq_record.annotations.items():
        if key == "comment" and not value:
            continue
        if isinstance(value, list):
            attributes[key] = value
        else:
            attributes[key] = [value]
    return attributes


def _strand(feature: object) -> str:
    strand = feature.location.strand
    if strand is None:
        return "."
    if strand < 0:
        return "-"
    return "+"


def _escape(value: object) -> str:
    return quote(str(value).rstrip(), safe=": /")
