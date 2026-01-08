#!/usr/bin/env python3
"""
Convert GenBank files to JSON format with comprehensive metadata and features.

This module reads GenBank files (containing one or multiple records/contigs) using Biopython
and writes a single JSON file matching a specific schema compatible with Bakta JSON format.

Coordinate System:
    Uses 1-based inclusive coordinates (start, stop) for features.
    BioPython's 0-based half-open [start, end) is converted to 1-based [start, stop].

Frame Convention:
    Frame is 0-based (0, 1, 2) derived from GenBank's codon_start (1, 2, 3).

GC Content:
    Calculated as (G+C) / (A+C+G+T), ignoring N and other ambiguous bases.

Hash Digest:
    MD5 hexdigest is used for amino acid sequences (Bakta compatible).
"""
import logging
import sys
import json
import argparse
import hashlib
from typing import List, Dict, Any, Optional
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

try:
    from Bio.SeqUtils import molecular_weight
    from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
    HAS_SEQ_UTILS = True
except ImportError:
    HAS_SEQ_UTILS = False

try:
    from .bacteria import gram_positive, gram_negative
except ImportError:
    # If running as standalone script
    try:
        from bacteria import gram_positive, gram_negative
    except ImportError:
        # Define minimal sets if bacteria.py not available
        gram_positive = set()
        gram_negative = set()

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


def calculate_n50(lengths: List[int]) -> int:
    """
    Calculate N50 from a list of contig lengths.
    
    N50 is the sequence length of the shortest contig at 50% of the total genome length.
    
    Args:
        lengths: List of contig lengths
        
    Returns:
        N50 value
    """
    if not lengths:
        return 0
    
    sorted_lengths = sorted(lengths, reverse=True)
    total_length = sum(sorted_lengths)
    half_length = total_length / 2
    
    cumulative_length = 0
    for length in sorted_lengths:
        cumulative_length += length
        if cumulative_length >= half_length:
            return length
    
    return sorted_lengths[-1] if sorted_lengths else 0


def calculate_gc_content(sequences: List[str]) -> float:
    """
    Calculate GC content across all sequences.
    
    GC fraction is calculated as (G+C) / (A+C+G+T).
    N bases and other ambiguous bases are ignored in the denominator.
    
    Args:
        sequences: List of DNA sequences
        
    Returns:
        GC fraction (0.0 to 1.0)
    """
    gc_count = 0
    total_count = 0
    
    for seq in sequences:
        seq_upper = seq.upper()
        gc_count += seq_upper.count('G') + seq_upper.count('C')
        total_count += (seq_upper.count('A') + seq_upper.count('C') + 
                       seq_upper.count('G') + seq_upper.count('T'))
    
    return gc_count / total_count if total_count > 0 else 0.0


def calculate_n_ratio(sequences: List[str]) -> float:
    """
    Calculate the fraction of N bases across all sequences.
    
    Args:
        sequences: List of DNA sequences
        
    Returns:
        N ratio (0.0 to 1.0)
    """
    n_count = 0
    total_count = 0
    
    for seq in sequences:
        seq_upper = seq.upper()
        n_count += seq_upper.count('N')
        total_count += len(seq)
    
    return n_count / total_count if total_count > 0 else 0.0


def extract_genome_metadata(records: List[SeqRecord], args: argparse.Namespace) -> Dict[str, Any]:
    """
    Extract genome-level metadata from records and command-line arguments.
    
    Priority: command-line args > GenBank annotations > gram stain inference
    
    Args:
        records: List of SeqRecord objects
        args: Command-line arguments
        
    Returns:
        Genome metadata dictionary
    """
    genome = {
        "genus": "NA",
        "species": "NA",
        "strain": "NA",
        "complete": False,
        "gram": "?",
        "translation_table": 11  # Default bacterial translation table
    }
    
    # Try to extract from first record annotations
    if records:
        record = records[0]
        annotations = record.annotations
        
        # Extract organism information
        if 'organism' in annotations:
            organism = annotations['organism']
            parts = organism.split()
            if len(parts) >= 1:
                genome['genus'] = parts[0]
            if len(parts) >= 2:
                genome['species'] = parts[1]
            if len(parts) >= 3:
                genome['strain'] = ' '.join(parts[2:])
        
        # Check for completeness
        if 'comment' in annotations:
            comment = annotations['comment'].lower()
            if 'complete' in comment:
                genome['complete'] = True
        
        # Check description for completeness
        if record.description and 'complete genome' in record.description.lower():
            genome['complete'] = True
        
        # Check topology for completeness indication
        for rec in records:
            if rec.annotations.get('topology', '').lower() == 'circular':
                genome['complete'] = True
                break
        
        # Infer gram stain from genus if not provided
        if genome['genus'] != "NA" and not args.gram:
            if genome['genus'] in gram_positive:
                genome['gram'] = '+'
            elif genome['genus'] in gram_negative:
                genome['gram'] = '-'
    
    # Override with command-line arguments if provided
    if args.genus:
        genome['genus'] = args.genus
    if args.species:
        genome['species'] = args.species
    if args.strain:
        genome['strain'] = args.strain
    if args.gram:
        genome['gram'] = args.gram
    if args.translation_table:
        genome['translation_table'] = args.translation_table
    
    return genome


def create_sequence_object(record: SeqRecord) -> Dict[str, Any]:
    """
    Create a sequence object from a GenBank record.
    
    Coordinates are 1-based inclusive.
    
    Args:
        record: BioPython SeqRecord
        
    Returns:
        Sequence object dictionary
    """
    seq_obj = {
        "id": record.id,
        "description": record.description,
        "sequence": str(record.seq).upper(),
        "length": len(record.seq),
        "complete": False,
        "type": "sequence",
        "topology": "linear",
        "simple_id": record.id.split('.')[0]  # Strip version number if present
    }
    
    # Determine topology
    if 'topology' in record.annotations:
        seq_obj['topology'] = record.annotations['topology']
    
    # Determine completeness
    if 'comment' in record.annotations:
        comment = record.annotations['comment'].lower()
        if 'complete' in comment:
            seq_obj['complete'] = True
    
    if seq_obj['topology'].lower() == 'circular':
        seq_obj['complete'] = True
    
    # Determine type (chromosome, plasmid, or contig)
    description_lower = record.description.lower()
    if 'chromosome' in description_lower:
        seq_obj['type'] = 'chromosome'
    elif 'plasmid' in description_lower:
        seq_obj['type'] = 'plasmid'
    elif 'contig' in description_lower:
        seq_obj['type'] = 'contig'
    
    # Add name if it's meaningful and different from id
    if record.name and record.name != record.id and record.name != '<unknown name>':
        seq_obj['name'] = record.name
    
    return seq_obj


def extract_feature_sequence(feature: SeqFeature, record: SeqRecord) -> str:
    """
    Extract nucleotide sequence for a feature (handles compound locations).
    
    Args:
        feature: BioPython SeqFeature
        record: BioPython SeqRecord
        
    Returns:
        Nucleotide sequence string
    """
    try:
        return str(feature.extract(record.seq))
    except (ValueError, AttributeError) as e:
        # Log or handle specific extraction errors
        return ""


def translate_sequence(nt_seq: str, translation_table: int) -> str:
    """
    Translate nucleotide sequence to amino acids.
    
    Args:
        nt_seq: Nucleotide sequence
        translation_table: NCBI translation table number
        
    Returns:
        Amino acid sequence
    """
    try:
        # Remove incomplete codons at the end
        complete_codons = (len(nt_seq) // 3) * 3
        if complete_codons > 0:
            seq = Seq(nt_seq[:complete_codons])
            aa_seq = str(seq.translate(table=translation_table))
            # Remove trailing stop codons
            return aa_seq.rstrip('*')
        return ""
    except (ValueError, KeyError) as e:
        # Handle invalid translation table or invalid nucleotide sequence
        return ""


def create_feature_object(feature: SeqFeature, record: SeqRecord, 
                          translation_table: int) -> Optional[Dict[str, Any]]:
    """
    Create a feature object from a GenBank feature.
    
    Coordinates are 1-based inclusive (start is +1, stop is unchanged from BioPython).
    
    Args:
        feature: BioPython SeqFeature
        record: BioPython SeqRecord
        translation_table: Translation table for CDS
        
    Returns:
        Feature object dictionary or None if not a relevant feature type
    """
    # Only process relevant feature types
    relevant_types = ['CDS', 'tRNA', 'rRNA', 'tmRNA', 'ncRNA']
    if feature.type not in relevant_types:
        return None
    
    qualifiers = feature.qualifiers
    
    # Determine strand
    if feature.location.strand is None:
        strand = 0
    elif feature.location.strand >= 0:
        strand = 1
    else:
        strand = -1
    
    # Base feature object
    # Coordinates: BioPython uses 0-based half-open intervals [start, end)
    # We convert to 1-based inclusive: start+1, end (end is already correct)
    feat_obj = {
        "type": feature.type,
        "contig": record.id,
        "start": int(feature.location.start) + 1,  # Convert to 1-based
        "stop": int(feature.location.end),  # Already correct for 1-based inclusive
        "strand": strand
    }
    
    # Extract common qualifiers
    if 'product' in qualifiers:
        feat_obj['product'] = qualifiers['product'][0]
    
    if 'gene' in qualifiers:
        feat_obj['gene'] = qualifiers['gene'][0]
    
    # Extract locus (try multiple qualifier names)
    locus = None
    for locus_key in ['locus_tag', 'protein_id', 'old_locus_tag']:
        if locus_key in qualifiers:
            locus = qualifiers[locus_key][0]
            break
    if locus:
        feat_obj['locus'] = locus
    
    # Extract db_xrefs and sort them
    if 'db_xref' in qualifiers:
        feat_obj['db_xrefs'] = sorted(qualifiers['db_xref'])
    
    # Extract nucleotide sequence
    nt_seq = extract_feature_sequence(feature, record)
    if nt_seq:
        feat_obj['nt'] = nt_seq
    
    # CDS-specific handling
    if feature.type == 'CDS':
        # Extract frame from codon_start (1,2,3 -> 0,1,2)
        if 'codon_start' in qualifiers:
            try:
                codon_start = int(qualifiers['codon_start'][0])
                feat_obj['frame'] = codon_start - 1  # Convert 1-based to 0-based
            except (ValueError, IndexError):
                # If codon_start is missing or malformed, fall back to the default frame (0)
                pass
        
        # Extract or translate amino acid sequence
        if 'translation' in qualifiers:
            aa_seq = qualifiers['translation'][0]
        elif nt_seq:
            # Apply frame offset if present
            frame_offset = feat_obj.get('frame', 0)
            aa_seq = translate_sequence(nt_seq[frame_offset:], translation_table)
        else:
            aa_seq = ""
        
        if aa_seq:
            feat_obj['aa'] = aa_seq
            # Calculate MD5 hexdigest (Bakta compatible)
            feat_obj['aa_hexdigest'] = hashlib.md5(aa_seq.encode()).hexdigest()
            
            # Calculate protein sequence statistics
            if HAS_SEQ_UTILS:
                try:
                    seq_stats = {}
                    seq_stats['molecular_weight'] = molecular_weight(aa_seq, seq_type="protein")
                    seq_stats['isoelectric_point'] = IsoelectricPoint(aa_seq).pi()
                    feat_obj['seq_stats'] = seq_stats
                except Exception:
                    # SeqUtils failures should not break JSON export; omit optional seq_stats on error.
                    pass
        
        # Check if hypothetical
        product = feat_obj.get('product', '')
        if 'hypothetical' in product.lower():
            feat_obj['hypothetical'] = True
        
        # Extract start_type based on actual start codon, if available
        start_type = "NA"
        if nt_seq:
            # Use frame offset (0, 1, 2) to determine the first codon
            frame_offset = feat_obj.get('frame', 0)
            try:
                frame_offset_int = int(frame_offset)
            except (TypeError, ValueError):
                frame_offset_int = 0
            if len(nt_seq) >= frame_offset_int + 3:
                start_codon = nt_seq[frame_offset_int:frame_offset_int + 3]
                start_type = str(start_codon).upper()
        feat_obj['start_type'] = start_type
        
        # Extract rbs_motif
        if 'ribosome_binding_site' in qualifiers:
            feat_obj['rbs_motif'] = qualifiers['ribosome_binding_site'][0]
        else:
            feat_obj['rbs_motif'] = "NA"
        
        # Better truncation detection
        # Use an unstripped translation (including stop codon) when possible
        raw_aa_seq: Optional[str] = None
        if nt_seq:
            frame_offset = feat_obj.get('frame', 0)
            try:
                raw_aa_seq = str(Seq(nt_seq[frame_offset:]).translate(table=translation_table))
            except Exception:
                # If translation fails for any reason, fall back to aa_seq-based heuristics only
                raw_aa_seq = None
        truncated_parts = []
        if aa_seq:
            if not aa_seq.startswith('M'):
                truncated_parts.append('5-prime')
            # 3-prime truncation detection:
            # - Prefer the raw translation (including stop codon) if available.
            #   A missing terminal '*' indicates a 3-prime truncation.
            # - If raw translation is unavailable, fall back to checking aa_seq
            #   (this may still contain '*' when taken directly from qualifiers).
            if raw_aa_seq is not None:
                if not raw_aa_seq.endswith('*'):
                    truncated_parts.append('3-prime')
            else:
                if aa_seq.endswith('*'):
                    truncated_parts.append('3-prime')
            
            if truncated_parts:
                # Join multiple truncation types if both exist
                feat_obj['truncated'] = '-'.join(truncated_parts)
        
        # Check for truncated/pseudo in qualifiers
        if ('truncated' in qualifiers or 'pseudo' in qualifiers) and 'truncated' not in feat_obj:
            feat_obj['truncated'] = 'unknown'

    try:
        # tRNA-specific handling
        if feature.type == 'tRNA':
            if 'product' in qualifiers:
                # Try to extract amino acid from product
                product = qualifiers['product'][0]
                if 'tRNA-' in product:
                    aa = product.split('tRNA-')[1].split()[0]
                    feat_obj['amino_acid'] = aa

            # Extract anticodon information
            if 'anticodon' in qualifiers:
                anticodon_info = qualifiers['anticodon'][0]
                # Format is typically like "(pos:1234..1236,aa:Met,seq:cat)"
                if 'seq:' in anticodon_info:
                    seq_part = anticodon_info.split('seq:')[1].rstrip(')')
                    feat_obj['anti_codon'] = seq_part
                if 'pos:' in anticodon_info:
                    pos_part = anticodon_info.split('pos:')[1].split(',')[0]
                    # Extract start position
                    if '..' in pos_part:
                        pos = int(pos_part.split('..')[0])
                        feat_obj['anti_codon_pos'] = pos
    except (IndexError, ValueError):
        # Malformed tRNA product/anticodon qualifiers are ignored
        logging.warning(f"Feature '{feature}' has a malformed tRNA or anti-codon product.")
        pass

    # Optional fields that might be in qualifiers
    if 'label' in qualifiers:
        feat_obj['label'] = qualifiers['label'][0]
    
    if 'note' in qualifiers:
        note = qualifiers['note'][0]
        # Try to extract class, score, evalue from notes if present
        if 'class:' in note.lower():
            feat_obj['class'] = note
    
    return feat_obj


def calculate_coding_ratio_accurate(records: List[SeqRecord]) -> float:
    """
    Calculate coding ratio by marking all coding bases in the genome.
    
    This method creates a mask for each base and marks it as coding if it's
    part of any CDS feature, providing a more accurate coding density measure.
    
    Args:
        records: List of SeqRecord objects
        
    Returns:
        Coding ratio (0.0 to 1.0)
    """
    total_coding = 0
    total_bases = 0
    
    for record in records:
        seq_len = len(record.seq)
        total_bases += seq_len
        coding_mask = [False] * seq_len
        
        for feature in record.features:
            if feature.type != 'CDS':
                continue
            
            start = int(feature.location.start)  # 0-based inclusive
            end = int(feature.location.end)      # 0-based exclusive
            
            # Mark all bases covered by this CDS as coding
            for i in range(start, end):
                if i < seq_len:  # Safety check
                    coding_mask[i] = True
        
        total_coding += sum(coding_mask)
    
    return total_coding / total_bases if total_bases > 0 else 0.0


def calculate_stats(records: List[SeqRecord], features: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Calculate genome statistics.
    
    Args:
        records: List of SeqRecord objects
        features: List of feature dictionaries
        
    Returns:
        Stats dictionary
    """
    sequences = [str(record.seq).upper() for record in records]
    lengths = [len(record.seq) for record in records]
    
    total_size = sum(lengths)
    
    # Calculate coding ratio using accurate base-by-base method
    coding_ratio = calculate_coding_ratio_accurate(records)
    
    stats = {
        "no_sequences": len(records),
        "size": total_size,
        "gc": calculate_gc_content(sequences),
        "n_ratio": calculate_n_ratio(sequences),
        "n50": calculate_n50(lengths),
        "coding_ratio": coding_ratio
    }
    
    return stats


def validate_json_output(data: Dict[str, Any]) -> None:
    """
    Validate the JSON output structure.
    
    Args:
        data: The complete JSON data structure
        
    Raises:
        AssertionError: If validation fails
    """
    # Check top-level keys
    required_keys = ['genome', 'stats', 'features', 'sequences', 'version']
    for key in required_keys:
        assert key in data, f"Missing top-level key: {key}"
    
    # Check sequences length matches stats
    assert len(data['sequences']) == data['stats']['no_sequences'], \
        "Number of sequences doesn't match no_sequences stat"
    
    # Check that every feature contig exists in sequences
    sequence_ids = {seq['id'] for seq in data['sequences']}
    for i, feat in enumerate(data['features']):
        assert 'contig' in feat, f"Feature {i} missing contig field"
        assert feat['contig'] in sequence_ids, \
            f"Feature {i} has contig '{feat['contig']}' not in sequences"


def genbank_to_json(genbank_path: str, args: argparse.Namespace) -> Dict[str, Any]:
    """
    Convert a GenBank file to JSON format.
    
    Args:
        genbank_path: Path to GenBank file
        args: Command-line arguments
        
    Returns:
        Complete JSON data structure
    """
    # Parse all records
    records = list(SeqIO.parse(genbank_path, "genbank"))
    
    if not records:
        raise ValueError(f"No records found in {genbank_path}")
    
    # Extract genome metadata
    genome = extract_genome_metadata(records, args)
    
    # Extract sequences
    sequences = [create_sequence_object(record) for record in records]
    
    # Extract features from all records
    features = []
    for record in records:
        for feature in record.features:
            feat_obj = create_feature_object(feature, record, genome['translation_table'])
            if feat_obj:
                features.append(feat_obj)
    
    # Sort features by contig, then start, then stop for deterministic output
    features.sort(key=lambda f: (f['contig'], f['start'], f['stop']))
    
    # Calculate stats
    stats = calculate_stats(records, features)
    
    # Version information
    version = {
        "bakta": args.bakta_version if args.bakta_version else "NA",
        "db": args.db_version if args.db_version else "NA"
    }
    
    # Assemble final JSON structure
    json_data = {
        "genome": genome,
        "stats": stats,
        "features": features,
        "sequences": sequences,
        "version": version
    }
    
    # Validate output
    validate_json_output(json_data)
    
    return json_data


def main():
    """
    Main entry point for command-line tool.
    """
    parser = argparse.ArgumentParser(
        description="Convert GenBank files to JSON format with comprehensive metadata"
    )
    
    # Required arguments
    parser.add_argument(
        '--genbank',
        required=True,
        help='Path to input GenBank file (.gbk or .gbff)'
    )
    parser.add_argument(
        '--out',
        required=True,
        help='Path to output JSON file'
    )
    
    # Optional version arguments
    parser.add_argument(
        '--bakta-version',
        default=None,
        help='Bakta version string (default: NA)'
    )
    parser.add_argument(
        '--db-version',
        default=None,
        help='Database version string (default: NA)'
    )
    
    # Optional genome metadata arguments
    parser.add_argument(
        '--genus',
        default=None,
        help='Genus name (overrides GenBank annotation)'
    )
    parser.add_argument(
        '--species',
        default=None,
        help='Species name (overrides GenBank annotation)'
    )
    parser.add_argument(
        '--strain',
        default=None,
        help='Strain name (overrides GenBank annotation)'
    )
    parser.add_argument(
        '--gram',
        default=None,
        choices=['+', '-'],
        help='Gram stain (+ or -)'
    )
    parser.add_argument(
        '--translation-table',
        type=int,
        default=None,
        help='NCBI translation table number (default: 11)'
    )
    
    args = parser.parse_args()
    
    try:
        # Convert GenBank to JSON
        json_data = genbank_to_json(args.genbank, args)
        
        # Write JSON output
        with open(args.out, 'w') as f:
            json.dump(json_data, f, indent=2)
        
        print(f"Successfully converted {args.genbank} to {args.out}")
        print(f"  Sequences: {json_data['stats']['no_sequences']}")
        print(f"  Features: {len(json_data['features'])}")
        print(f"  Total size: {json_data['stats']['size']} bp")
        
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)


if __name__ == '__main__':
    # Self-test mode if no arguments provided
    if len(sys.argv) == 1:
        print("GenBank to JSON Converter")
        print("=" * 50)
        print("\nUsage:")
        print("  python genbank_to_json.py --genbank input.gbk --out output.json")
        print("\nFor self-test, provide a GenBank file:")
        print("  python genbank_to_json.py --genbank test.gbk --out test.json")
        sys.exit(0)
    
    main()
