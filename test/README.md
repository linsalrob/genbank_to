# Test Suite for genbank_to

This directory contains the test suite for the genbank_to command-line tool.

## Test Files

- `NC_001417.gbk` - Test input file (GenBank format)
- `outputs/` - Directory containing expected output files for all format conversions
- `test_genbank_to.py` - Test suite that verifies all conversion formats

## Running the Tests

### Install Test Dependencies

```bash
pip install pytest
```

Or install with test extras:

```bash
pip install -e ".[test]"
```

### Run All Tests

```bash
pytest test/test_genbank_to.py -v
```

### Run Specific Tests

```bash
# Test PTT conversion only
pytest test/test_genbank_to.py::TestGenbankTo::test_ptt_conversion -v

# Test AMR conversion only
pytest test/test_genbank_to.py::TestGenbankTo::test_amr_conversion -v
```

## What the Tests Do

The test suite runs the following commands and verifies that the outputs match the expected reference files:

```bash
genbank_to -g test/NC_001417.gbk --ptt test/outputs/NC_001417.ptt
genbank_to -g test/NC_001417.gbk -a test/outputs/NC_001417.faa
genbank_to -g test/NC_001417.gbk -n test/outputs/NC_001417.fna
genbank_to -g test/NC_001417.gbk --orfs test/outputs/NC_001417.orfs
genbank_to -g test/NC_001417.gbk --functions test/outputs/NC_001417.functions
genbank_to -g test/NC_001417.gbk --amr test/outputs/NC_001417.amr
genbank_to -g test/NC_001417.gbk --gff3 test/outputs/NC_001417.gff3
genbank_to -g test/NC_001417.gbk --phage_finder test/outputs/NC_001417.phage_finder
genbank_to -g test/NC_001417.gbk --bakta-json test/outputs/NC_001417.json
```

Each test:
1. Runs the genbank_to command with the appropriate options
2. Generates output in a temporary directory
3. Compares the generated output with the expected reference file
4. Fails if the outputs don't match

## Test Coverage

The test suite covers all major output formats:
- PTT (protein table)
- FAA (amino acid sequences)
- FNA (nucleotide sequences)
- ORFs (open reading frames)
- Functions (protein ID and function table)
- AMR (AMRFinderPlus format - creates .gff, .faa, and .fna files)
- GFF3 (General Feature Format version 3)
- Phage Finder format
- Bakta JSON format
