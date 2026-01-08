Welcome to genbank_to's documentation!
=======================================

.. image:: https://img.shields.io/badge/Bioinformatics-EdwardsLab-03A9F4
   :target: https://edwards.flinders.edu.au/
   :alt: Edwards Lab

.. image:: https://www.zenodo.org/badge/481464683.svg
   :target: https://www.zenodo.org/badge/latestdoi/481464683
   :alt: DOI

.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License: MIT

.. image:: https://img.shields.io/pypi/pyversions/genbank_to.svg?style=flat-square&label=PyPi%20Versions
   :target: https://pypi.org/project/genbank_to/
   :alt: PyPi

**genbank_to** is a straightforward Python application and library for converting NCBI GenBank format files to a variety of other formats commonly used in bioinformatics workflows.

Overview
--------

The tool reads NCBI GenBank format files and converts them to various output formats including:

- FASTA nucleotide sequences (genome, ORFs)
- FASTA amino acid sequences (proteins)
- GFF3 format
- NCBI PTT format
- Function tables
- Bakta JSON format
- AMRFinderPlus format
- Phage Finder format

Both command-line and library interfaces are provided, making it easy to integrate into scripts and pipelines.

Key Features
------------

- **Multiple output formats**: Convert to numerous formats in a single command
- **Flexible input**: Handles single and multi-record GenBank files
- **Python library**: Import and use functions in your own scripts
- **Well-tested**: Comprehensive test suite ensures reliability
- **Active development**: Maintained by the Edwards Lab at Flinders University

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   quickstart
   usage
   output_formats
   examples
   api
   contributing
   changelog

Quick Example
-------------

Convert a GenBank file to multiple formats:

.. code-block:: bash

   genbank_to -g genome.gbk \
       -n genome.fna \
       -a proteins.faa \
       -o orfs.fna \
       --gff3 genome.gff3

Use as a Python library:

.. code-block:: python

   from GenBankToLib import genbank_to_faa, genbank_to_fna
   
   # Extract protein sequences
   for seqid, protid, sequence in genbank_to_faa('genome.gbk'):
       print(f">{protid}\n{sequence}")

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
