Quick Start
===========

This guide will get you up and running with genbank_to quickly.

Basic Usage
-----------

The most common use case is converting a GenBank file to FASTA format:

.. code-block:: bash

   genbank_to -g input.gbk -n output.fna

This reads ``input.gbk`` and writes the nucleotide sequence to ``output.fna``.

Common Conversions
------------------

Extract Genome Sequence
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   genbank_to -g genome.gbk -n genome.fna

Extract Protein Sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   genbank_to -g genome.gbk -a proteins.faa

Extract ORF Sequences (DNA)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   genbank_to -g genome.gbk -o orfs.fna

Generate GFF3 File
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   genbank_to -g genome.gbk --gff3 genome.gff3

Multiple Outputs at Once
~~~~~~~~~~~~~~~~~~~~~~~~~

You can request multiple output formats in a single command:

.. code-block:: bash

   genbank_to -g genome.gbk \
       -n genome.fna \
       -a proteins.faa \
       -o orfs.fna \
       --gff3 genome.gff3

Using as a Python Library
--------------------------

You can also use genbank_to as a library in your Python scripts:

Extract All Protein Sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from GenBankToLib import genbank_to_faa
   
   for seqid, protid, sequence in genbank_to_faa('genome.gbk'):
       print(f">{protid}")
       print(sequence)

Extract Genome Sequences
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from GenBankToLib import genbank_to_fna
   
   for seqid, sequence in genbank_to_fna('genome.gbk'):
       print(f">{seqid}")
       print(sequence)

Extract Functions
~~~~~~~~~~~~~~~~~

.. code-block:: python

   from GenBankToLib import genbank_to_functions
   
   for protid, function in genbank_to_functions('genome.gbk'):
       print(f"{protid}\t{function}")

Convert to JSON Format
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from GenBankToLib import genbank_to_json
   import json
   
   genome_info = {'gram': None, 'translation_table': 11}
   json_data = genbank_to_json('genome.gbk', genome_info)
   
   # Save to file
   with open('genome.json', 'w') as f:
       json.dump(json_data, f, indent=2)

Example Workflow
----------------

Here's a complete example workflow for analyzing a phage genome:

1. **Download a phage genome from NCBI**:

.. code-block:: bash

   # Example: Enterobacteria phage phiX174
   # (This is the test file in the repository)
   wget https://raw.githubusercontent.com/linsalrob/genbank_to/main/test/NC_001417.gbk

2. **Extract multiple formats**:

.. code-block:: bash

   genbank_to -g NC_001417.gbk \
       -n NC_001417.fna \
       -a NC_001417.faa \
       -o NC_001417_orfs.fna \
       -f NC_001417_functions.tsv \
       --gff3 NC_001417.gff3

3. **View the results**:

.. code-block:: bash

   # View genome sequence
   head NC_001417.fna
   
   # Count proteins
   grep -c ">" NC_001417.faa
   
   # View functions
   head NC_001417_functions.tsv

Working with Multi-GenBank Files
---------------------------------

If your GenBank file contains multiple sequences (separated by ``//``), you can split them:

.. code-block:: bash

   genbank_to -g multi.gbk --separate -n output

This creates separate files: ``output.seqid1.fna``, ``output.seqid2.fna``, etc.

Filtering by Sequence ID
~~~~~~~~~~~~~~~~~~~~~~~~~

Extract only specific sequences:

.. code-block:: bash

   genbank_to -g multi.gbk -i NC_001417 -n output.fna

Next Steps
----------

- Read the :doc:`usage` guide for detailed command-line options
- Explore the :doc:`output_formats` for all available formats
- Check the :doc:`examples` for more complex use cases
- Review the :doc:`api` documentation for library usage
