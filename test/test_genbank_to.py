"""
Test suite for genbank_to command-line tool.

This test suite runs genbank_to commands and verifies that the outputs
match the expected reference files in test/outputs/.
"""

import os
import subprocess
import tempfile
import filecmp
from pathlib import Path

import pytest


# Get the base directory for tests
TEST_DIR = Path(__file__).parent
INPUT_FILE = TEST_DIR / "NC_001417.gbk"
EXPECTED_OUTPUT_DIR = TEST_DIR / "outputs"


class TestGenbankTo:
    """Test suite for genbank_to conversions."""

    def run_genbank_to(self, *args):
        """Helper method to run genbank_to command."""
        cmd = ["genbank_to"] + list(args)
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False
        )
        if result.returncode != 0:
            pytest.fail(
                f"Command failed: {' '.join(cmd)}\n"
                f"stdout: {result.stdout}\n"
                f"stderr: {result.stderr}"
            )
        return result

    def compare_files(self, expected_file, actual_file):
        """Compare two files and assert they are identical."""
        assert os.path.exists(expected_file), f"Expected file not found: {expected_file}"
        assert os.path.exists(actual_file), f"Actual file not found: {actual_file}"
        
        # Use filecmp for binary comparison
        if not filecmp.cmp(expected_file, actual_file, shallow=False):
            # If files differ, show the difference
            with open(expected_file, 'r') as f1, open(actual_file, 'r') as f2:
                expected_content = f1.read()
                actual_content = f2.read()
            pytest.fail(
                f"Files differ:\n"
                f"Expected: {expected_file}\n"
                f"Actual: {actual_file}\n"
                f"Expected content:\n{expected_content[:500]}\n...\n"
                f"Actual content:\n{actual_content[:500]}\n..."
            )

    def test_ptt_conversion(self, tmp_path):
        """Test conversion to PTT format."""
        output_file = tmp_path / "NC_001417.ptt"
        expected_file = EXPECTED_OUTPUT_DIR / "NC_001417.ptt"
        
        self.run_genbank_to("-g", str(INPUT_FILE), "--ptt", str(output_file))
        self.compare_files(expected_file, output_file)

    def test_faa_conversion(self, tmp_path):
        """Test conversion to FAA (amino acid) format."""
        output_file = tmp_path / "NC_001417.faa"
        expected_file = EXPECTED_OUTPUT_DIR / "NC_001417.faa"
        
        self.run_genbank_to("-g", str(INPUT_FILE), "-a", str(output_file))
        self.compare_files(expected_file, output_file)

    def test_fna_conversion(self, tmp_path):
        """Test conversion to FNA (nucleotide) format."""
        output_file = tmp_path / "NC_001417.fna"
        expected_file = EXPECTED_OUTPUT_DIR / "NC_001417.fna"
        
        self.run_genbank_to("-g", str(INPUT_FILE), "-n", str(output_file))
        self.compare_files(expected_file, output_file)

    def test_orfs_conversion(self, tmp_path):
        """Test conversion to ORFs format."""
        output_file = tmp_path / "NC_001417.orfs"
        expected_file = EXPECTED_OUTPUT_DIR / "NC_001417.orfs"
        
        self.run_genbank_to("-g", str(INPUT_FILE), "--orfs", str(output_file))
        self.compare_files(expected_file, output_file)

    def test_functions_conversion(self, tmp_path):
        """Test conversion to functions format."""
        output_file = tmp_path / "NC_001417.functions"
        expected_file = EXPECTED_OUTPUT_DIR / "NC_001417.functions"
        
        self.run_genbank_to("-g", str(INPUT_FILE), "--functions", str(output_file))
        self.compare_files(expected_file, output_file)

    def test_amr_conversion(self, tmp_path):
        """Test conversion to AMR format (creates multiple files)."""
        # AMR format creates three files: .amr.gff, .amr.faa, .amr.fna
        base_output = tmp_path / "NC_001417.amr"
        
        self.run_genbank_to("-g", str(INPUT_FILE), "--amr", str(base_output))
        
        # Check all three AMR output files
        expected_files = [
            "NC_001417.amr.gff",
            "NC_001417.amr.faa",
            "NC_001417.amr.fna"
        ]
        
        for filename in expected_files:
            expected_file = EXPECTED_OUTPUT_DIR / filename
            actual_file = tmp_path / filename
            self.compare_files(expected_file, actual_file)

    def test_gff3_conversion(self, tmp_path):
        """Test conversion to GFF3 format."""
        output_file = tmp_path / "NC_001417.gff3"
        expected_file = EXPECTED_OUTPUT_DIR / "NC_001417.gff3"
        
        self.run_genbank_to("-g", str(INPUT_FILE), "--gff3", str(output_file))
        self.compare_files(expected_file, output_file)

    def test_phage_finder_conversion(self, tmp_path):
        """Test conversion to phage_finder format."""
        output_file = tmp_path / "NC_001417.phage_finder"
        expected_file = EXPECTED_OUTPUT_DIR / "NC_001417.phage_finder"
        
        self.run_genbank_to("-g", str(INPUT_FILE), "--phage_finder", str(output_file))
        self.compare_files(expected_file, output_file)

    def test_bakta_json_conversion(self, tmp_path):
        """Test conversion to Bakta JSON format."""
        output_file = tmp_path / "NC_001417.json"
        expected_file = EXPECTED_OUTPUT_DIR / "NC_001417.json"
        
        self.run_genbank_to("-g", str(INPUT_FILE), "--bakta-json", str(output_file))
        self.compare_files(expected_file, output_file)

    def test_input_file_exists(self):
        """Verify that the test input file exists."""
        assert INPUT_FILE.exists(), f"Test input file not found: {INPUT_FILE}"

    def test_expected_outputs_exist(self):
        """Verify that all expected output files exist."""
        expected_files = [
            "NC_001417.ptt",
            "NC_001417.faa",
            "NC_001417.fna",
            "NC_001417.orfs",
            "NC_001417.functions",
            "NC_001417.amr.gff",
            "NC_001417.amr.faa",
            "NC_001417.amr.fna",
            "NC_001417.gff3",
            "NC_001417.phage_finder",
            "NC_001417.json"
        ]
        
        for filename in expected_files:
            expected_file = EXPECTED_OUTPUT_DIR / filename
            assert expected_file.exists(), f"Expected output file not found: {expected_file}"
