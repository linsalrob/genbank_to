Changelog
=========

This page documents the version history and changes in genbank_to.

Version History
---------------

For the most up-to-date version information, see the `VERSION file <https://github.com/linsalrob/genbank_to/blob/main/VERSION>`_ in the repository.

Development Version
~~~~~~~~~~~~~~~~~~~

Current development happens on the `main branch <https://github.com/linsalrob/genbank_to>`_.

Previous Releases
~~~~~~~~~~~~~~~~~

For a complete list of releases and their changes, see the `GitHub Releases page <https://github.com/linsalrob/genbank_to/releases>`_.

Recent Updates
--------------

The genbank_to tool has been continuously improved with the following types of enhancements:

New Features
~~~~~~~~~~~~

- **Bakta JSON format**: Support for output compatible with Bakta genome annotation
- **AMRFinderPlus format**: Specialized output for antimicrobial resistance gene annotation
- **Phage Finder format**: Support for prophage identification workflows
- **GFF3 format**: Standard genomic feature format support
- **Multi-record support**: Handle GenBank files with multiple sequences
- **Sequence filtering**: Extract specific sequences by ID
- **Gram stain detection**: Automatic Gram stain determination for known genera
- **Complex headers**: Detailed FASTA headers with organism and annotation information

Improvements
~~~~~~~~~~~~

- **Memory efficiency**: Generator-based functions for processing large files
- **Error handling**: Better error messages and pseudogene handling
- **Gzip support**: Automatic detection and handling of compressed files
- **Logging**: Comprehensive logging for debugging
- **Type hints**: Better IDE support and code documentation
- **Testing**: Comprehensive test suite for reliability

Bug Fixes
~~~~~~~~~

- **BioPython compatibility**: Fixed issues with different BioPython versions
- **Pseudogene handling**: Proper handling of pseudogenes with optional inclusion
- **Coordinate systems**: Correct conversion between BioPython and GFF3 coordinates
- **Frame calculation**: Proper frame determination from codon_start
- **Translation handling**: Better handling of missing translations

Migration Notes
---------------

From NCBI Downloads to genbank_to
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you previously used NCBI's downloadable genome files (PTT, FAA, FNA), you can now generate these from GenBank files:

.. code-block:: bash

   # Old: Download from NCBI
   # wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/...
   
   # New: Convert from GenBank
   genbank_to -g genome.gbk \
       -n genome.fna \
       -a genome.faa \
       -p genome.ptt

From deprekate/genbank
~~~~~~~~~~~~~~~~~~~~~~~

If you're familiar with the `deprekate/genbank <https://github.com/deprekate/genbank>`_ package, genbank_to provides similar functionality with additional output formats:

.. code-block:: python

   # deprekate/genbank style
   import genbank
   gb = genbank.load('genome.gbk')
   
   # genbank_to style
   from GenBankToLib import genbank_to_faa, genbank_to_fna
   for seqid, sequence in genbank_to_fna('genome.gbk'):
       print(f">{seqid}\n{sequence}")

Citation
--------

If you use genbank_to in your research, please cite:

.. code-block:: bibtex

   @software{edwards_genbank_to,
     author       = {Edwards, Rob},
     title        = {genbank_to: Convert GenBank format files to other formats},
     year         = {2024},
     publisher    = {GitHub},
     url          = {https://github.com/linsalrob/genbank_to},
     doi          = {10.5281/zenodo.xxxxx}
   }

See the `CITATION.cff file <https://github.com/linsalrob/genbank_to/blob/main/CITATION.cff>`_ for the most current citation information.

Roadmap
-------

Planned Features
~~~~~~~~~~~~~~~~

We're always working to improve genbank_to. Planned features include:

- Additional output formats based on user requests
- Support for other input formats (GFF3, EMBL)
- Performance optimizations for very large genomes
- Parallel processing for batch conversions
- Integration with more bioinformatics pipelines

Request a Feature
~~~~~~~~~~~~~~~~~

Have an idea for a new feature? `Open an issue <https://github.com/linsalrob/genbank_to/issues/new>`_ on GitHub!

Contributing
~~~~~~~~~~~~

See the :doc:`contributing` guide for information on how to contribute to genbank_to.

Versioning
----------

genbank_to follows semantic versioning:

- **Major version** (X.0.0): Incompatible API changes
- **Minor version** (0.X.0): New features, backwards compatible
- **Patch version** (0.0.X): Bug fixes, backwards compatible

Support
-------

Getting Help
~~~~~~~~~~~~

If you need help:

1. Check the :doc:`usage` documentation
2. Look through :doc:`examples`
3. Search `existing issues <https://github.com/linsalrob/genbank_to/issues>`_
4. Create a new issue if needed

Reporting Issues
~~~~~~~~~~~~~~~~

Found a bug? Please report it with:

- genbank_to version (``genbank_to --version``)
- Python version
- Operating system
- Steps to reproduce
- Expected vs actual behavior

See :doc:`contributing` for more details.

Security
--------

If you discover a security vulnerability, please email raedwards@gmail.com directly rather than creating a public issue.

License
-------

genbank_to is released under the MIT License. See the `LICENSE file <https://github.com/linsalrob/genbank_to/blob/main/LICENSE>`_ for details.

Acknowledgments
---------------

genbank_to is developed by the Edwards Lab at Flinders University. We thank all contributors and users who provide feedback and suggestions.

Special thanks to:

- The BioPython team for their excellent library
- The NCBI for maintaining GenBank
- All our contributors on GitHub
- The bioinformatics community for feedback and feature requests
