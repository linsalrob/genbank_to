Command-Line Usage
==================

This page provides detailed documentation for all command-line options available in genbank_to.

Basic Syntax
------------

.. code-block:: bash

   genbank_to [OPTIONS]

Required Arguments
------------------

``-g``, ``--genbank`` FILENAME
   Path to the input GenBank file (required).

   .. code-block:: bash

      genbank_to -g genome.gbk -n output.fna

Output Format Options
---------------------

Nucleotide Outputs
~~~~~~~~~~~~~~~~~~

``-n``, ``--nucleotide`` FILENAME
   Output the complete nucleotide sequence(s) from the GenBank file (e.g., the genome sequence).

   .. code-block:: bash

      genbank_to -g genome.gbk -n genome.fna

``-o``, ``--orfs`` FILENAME
   Output the DNA sequences of all open reading frames (ORFs/CDS features).

   .. code-block:: bash

      genbank_to -g genome.gbk -o orfs.fna

Protein Outputs
~~~~~~~~~~~~~~~

``-a``, ``--aminoacids`` FILENAME
   Output the amino acid (protein) sequences for all coding sequences.

   .. code-block:: bash

      genbank_to -g genome.gbk -a proteins.faa

Complex Format Outputs
~~~~~~~~~~~~~~~~~~~~~~~

``-p``, ``--ptt`` FILENAME
   Output in NCBI PTT (Protein Table) format. This is a somewhat deprecated NCBI format from their genome downloads.

   .. code-block:: bash

      genbank_to -g genome.gbk -p genome.ptt

``-f``, ``--functions`` FILENAME
   Output a tab-separated table with protein ID and function (product) columns.

   .. code-block:: bash

      genbank_to -g genome.gbk -f functions.tsv

``--gff3`` FILENAME
   Output in GFF3 (Generic Feature Format version 3) format.

   .. code-block:: bash

      genbank_to -g genome.gbk --gff3 genome.gff3

``--amr`` BASENAME
   Output files in the format required by `NCBI AMRFinderPlus <https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus>`_.
   
   Creates three files:
   
   - ``BASENAME.gff`` - GFF format annotation
   - ``BASENAME.faa`` - Amino acid sequences
   - ``BASENAME.fna`` - Nucleotide sequences

   .. code-block:: bash

      genbank_to -g genome.gbk --amr genome_amr

``--phage_finder`` FILENAME
   Output in the format required by `phage_finder <http://phage-finder.sourceforge.net/>`_.

   .. code-block:: bash

      genbank_to -g phage.gbk --phage_finder phage.pf

``--bakta-json`` FILENAME
   Output JSON format similar to that created by `Bakta <https://github.com/oschwengers/bakta>`_.

   .. code-block:: bash

      genbank_to -g genome.gbk --bakta-json genome.json

Bakta JSON Metadata Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These options are only valid when ``--bakta-json`` is specified:

``--bakta-version`` STRING
   Bakta version string (default: NA). For recording which version of Bakta was used.

``--db-version`` STRING
   Database version string (default: NA). For recording the annotation database version.

``--genus`` STRING
   Genus name. Overrides the genus from GenBank annotation.

``--species`` STRING
   Species name. Overrides the species from GenBank annotation.

``--strain`` STRING
   Strain designation. Overrides the strain from GenBank annotation.

``--gram`` {+,-}
   Gram stain result (+ for Gram-positive, - for Gram-negative). If not provided, the tool will attempt to determine this from the genus name.

``--translation-table`` NUMBER
   NCBI translation table number (default: 11). Specify if you used a non-standard genetic code.

Example:

.. code-block:: bash

   genbank_to -g genome.gbk --bakta-json genome.json \
       --genus Escherichia \
       --species coli \
       --strain K-12 \
       --gram - \
       --translation-table 11

Output Modifiers
----------------

``-c``, ``--complex``
   Use complex/detailed identifier lines in the output. Includes additional information such as organism name, location, and product description in the FASTA headers.

   .. code-block:: bash

      genbank_to -g genome.gbk -a proteins.faa --complex

``--pseudo``
   Include pseudogenes in the output. By default, pseudogenes are skipped because they often cause BioPython errors. Use this flag to attempt including them.

   .. code-block:: bash

      genbank_to -g genome.gbk -a proteins.faa --pseudo

``-i``, ``--seqid`` ID
   Only output specific sequence ID(s). Can be specified multiple times to select multiple sequences. Automatically enables ``--separate``.

   .. code-block:: bash

      # Extract a single sequence
      genbank_to -g multi.gbk -i NC_001417 -n output.fna
      
      # Extract multiple sequences
      genbank_to -g multi.gbk -i NC_001417 -i NC_001418 -n output.fna

``--separate``
   Separate multi-record GenBank files into individual output files. Each sequence gets its own file with the sequence ID in the filename.

   .. code-block:: bash

      # Creates output.NC_001417.fna, output.NC_001418.fna, etc.
      genbank_to -g multi.gbk --separate -n output
      
      # With no other options, outputs separate GenBank files
      genbank_to -g multi.gbk --separate

``-z``, ``--zip``
   Compress the output using gzip. Experimental feature that may not work with all output formats.

   .. code-block:: bash

      genbank_to -g genome.gbk -f functions.tsv --zip

Logging and Debugging
---------------------

``--log`` FILENAME
   Specify the log file location (default: ``genbank_to.log``).

   .. code-block:: bash

      genbank_to -g genome.gbk -n output.fna --log my_log.txt

``-d``, ``--debug``
   Enable debug-level logging for troubleshooting.

   .. code-block:: bash

      genbank_to -g genome.gbk -n output.fna --debug

Other Options
-------------

``-v``, ``--version``
   Show the version number and exit.

   .. code-block:: bash

      genbank_to --version

``-h``, ``--help``
   Show help message and exit.

   .. code-block:: bash

      genbank_to --help

Complete Example
----------------

Here's a comprehensive example using multiple options:

.. code-block:: bash

   genbank_to \
       --genbank genome.gbk \
       --nucleotide genome.fna \
       --aminoacids proteins.faa \
       --orfs orfs.fna \
       --functions functions.tsv \
       --gff3 genome.gff3 \
       --bakta-json genome.json \
       --genus Escherichia \
       --species coli \
       --strain K-12 \
       --gram - \
       --complex \
       --log conversion.log \
       --debug

This command will:

1. Read the GenBank file ``genome.gbk``
2. Output the genome sequence to ``genome.fna``
3. Output protein sequences to ``proteins.faa`` with complex headers
4. Output ORF sequences to ``orfs.fna``
5. Output a function table to ``functions.tsv``
6. Output GFF3 format to ``genome.gff3``
7. Output Bakta JSON to ``genome.json`` with custom metadata
8. Write debug logs to ``conversion.log``
