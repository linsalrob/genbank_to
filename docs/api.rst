API Reference
=============

This page documents the Python library API for genbank_to. All functions can be imported from the ``GenBankToLib`` package.

.. code-block:: python

   from GenBankToLib import (
       genbank_to_faa,
       genbank_to_fna,
       genbank_to_orfs,
       genbank_to_functions,
       genbank_to_ptt,
       genbank_to_gff,
       genbank_to_phage_finder,
       genbank_to_amrfinder,
       genbank_seqio
   )

Core Functions
--------------

genbank_to_fna
~~~~~~~~~~~~~~

.. code-block:: python

   genbank_to_fna(gbkf, include_definition=False)

Extract nucleotide sequences from a GenBank file.

**Parameters:**

- ``gbkf`` (str): Path to the GenBank file.
- ``include_definition`` (bool, optional): If True, includes the GenBank definition line with the sequence ID. Default is False.

**Yields:**

- ``tuple``: (sequence_id, sequence) where:
  
  - ``sequence_id`` (str): The sequence identifier
  - ``sequence`` (Bio.Seq.Seq): The nucleotide sequence

**Example:**

.. code-block:: python

   from GenBankToLib import genbank_to_fna
   
   for seqid, sequence in genbank_to_fna('genome.gbk'):
       print(f">{seqid}")
       print(sequence)

genbank_to_faa
~~~~~~~~~~~~~~

.. code-block:: python

   genbank_to_faa(gbkf, complexheader=False, skip_pseudo=True)

Extract amino acid (protein) sequences from a GenBank file.

**Parameters:**

- ``gbkf`` (str): Path to the GenBank file.
- ``complexheader`` (bool, optional): If True, creates detailed headers with organism, location, product, and database references. Default is False.
- ``skip_pseudo`` (bool, optional): If True, skips pseudogenes. Default is True.

**Yields:**

- ``tuple``: (sequence_id, protein_id, sequence) where:
  
  - ``sequence_id`` (str): The parent sequence identifier
  - ``protein_id`` (str): The protein identifier (protein_id, locus_tag, or db_xref)
  - ``sequence`` (str): The amino acid sequence

**Example:**

.. code-block:: python

   from GenBankToLib import genbank_to_faa
   
   # Simple headers
   for seqid, protid, sequence in genbank_to_faa('genome.gbk'):
       print(f">{protid}")
       print(sequence)
   
   # Complex headers
   for seqid, protid, sequence in genbank_to_faa('genome.gbk', complexheader=True):
       print(f">{protid}")
       print(sequence)

genbank_to_orfs
~~~~~~~~~~~~~~~

.. code-block:: python

   genbank_to_orfs(gbkf, complexheader=False, skip_pseudo=True)

Extract DNA sequences of open reading frames (ORFs) from a GenBank file.

**Parameters:**

- ``gbkf`` (str): Path to the GenBank file.
- ``complexheader`` (bool, optional): If True, creates detailed headers. Default is False.
- ``skip_pseudo`` (bool, optional): If True, skips pseudogenes. Default is True.

**Yields:**

- ``tuple``: (sequence_id, protein_id, sequence) where:
  
  - ``sequence_id`` (str): The parent sequence identifier
  - ``protein_id`` (str): The protein/feature identifier
  - ``sequence`` (str): The DNA sequence of the ORF

**Example:**

.. code-block:: python

   from GenBankToLib import genbank_to_orfs
   
   for seqid, protid, sequence in genbank_to_orfs('genome.gbk'):
       print(f">{protid}")
       print(sequence)

genbank_to_functions
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   genbank_to_functions(gbkf, seqid=False, skip_pseudo=True)

Extract protein functions (products) from a GenBank file.

**Parameters:**

- ``gbkf`` (str): Path to the GenBank file.
- ``seqid`` (bool, optional): If True, includes the sequence ID in the output. Default is False.
- ``skip_pseudo`` (bool, optional): If True, skips pseudogenes. Default is True.

**Yields:**

- ``tuple``: If seqid=False: (protein_id, function)
- ``tuple``: If seqid=True: (sequence_id, protein_id, function)

**Example:**

.. code-block:: python

   from GenBankToLib import genbank_to_functions
   
   # Without sequence ID
   for protid, function in genbank_to_functions('genome.gbk'):
       print(f"{protid}\t{function}")
   
   # With sequence ID
   for seqid, protid, function in genbank_to_functions('genome.gbk', seqid=True):
       print(f"{seqid}\t{protid}\t{function}")

genbank_to_ptt
~~~~~~~~~~~~~~

.. code-block:: python

   genbank_to_ptt(gbkf, printout=False)

Convert GenBank file to NCBI Protein Table (PTT) format.

**Parameters:**

- ``gbkf`` (str): Path to the GenBank file.
- ``printout`` (bool, optional): If True, prints the table to stdout. Default is False.

**Returns:**

- ``list``: List of lists, where each inner list contains:
  [location, strand, length, gi, gene, synonym, cog, product]

**Example:**

.. code-block:: python

   from GenBankToLib import genbank_to_ptt
   
   table = genbank_to_ptt('genome.gbk')
   for row in table:
       print('\t'.join(map(str, row)))

genbank_to_gff
~~~~~~~~~~~~~~

.. code-block:: python

   genbank_to_gff(gbkf, out_gff)

Convert GenBank file to GFF3 format.

**Parameters:**

- ``gbkf`` (str): Path to the GenBank file.
- ``out_gff`` (str): Path to the output GFF3 file.

**Returns:**

- ``None``: Writes output directly to file.

**Example:**

.. code-block:: python

   from GenBankToLib import genbank_to_gff
   
   genbank_to_gff('genome.gbk', 'genome.gff3')

genbank_to_phage_finder
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   genbank_to_phage_finder(gbkf)

Convert GenBank file to phage_finder format.

**Parameters:**

- ``gbkf`` (str): Path to the GenBank file.

**Yields:**

- ``list``: [contig_id, contig_length, gene_id, start, end, function]

**Example:**

.. code-block:: python

   from GenBankToLib import genbank_to_phage_finder
   
   for row in genbank_to_phage_finder('genome.gbk'):
       print('\t'.join(map(str, row)))

genbank_to_amrfinder
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   genbank_to_amrfinder(gbkf, amrout)

Convert GenBank file to AMRFinderPlus format.

**Parameters:**

- ``gbkf`` (str): Path to the GenBank file.
- ``amrout`` (str): Base name for output files (creates .gff, .faa, and .fna files).

**Returns:**

- ``None``: Writes output directly to files.

**Example:**

.. code-block:: python

   from GenBankToLib import genbank_to_amrfinder
   
   genbank_to_amrfinder('genome.gbk', 'output')
   # Creates: output.gff, output.faa, output.fna

genbank_seqio
~~~~~~~~~~~~~

.. code-block:: python

   genbank_seqio(gbkf)

Get a BioPython SeqIO parser for the GenBank file.

**Parameters:**

- ``gbkf`` (str): Path to the GenBank file (can be gzipped).

**Returns:**

- ``tuple``: (parser, file_handle) where:
  
  - ``parser``: BioPython SeqIO parser object
  - ``file_handle``: File handle (should be closed when done)

**Example:**

.. code-block:: python

   from GenBankToLib import genbank_seqio
   
   parser, handle = genbank_seqio('genome.gbk')
   for record in parser:
       print(f"Processing {record.id}")
       # Do something with record
   handle.close()

Utility Functions
-----------------

feature_id
~~~~~~~~~~

.. code-block:: python

   feature_id(seq, feat)

Choose the appropriate identifier for a feature.

**Parameters:**

- ``seq`` (Bio.SeqRecord.SeqRecord): The sequence record.
- ``feat`` (Bio.SeqFeature.SeqFeature): The feature.

**Returns:**

- ``str``: The feature ID (protein_id, locus_tag, db_xref, or location-based ID).

**Example:**

.. code-block:: python

   from Bio import SeqIO
   from GenBankToLib.genbank import feature_id
   
   for record in SeqIO.parse('genome.gbk', 'genbank'):
       for feature in record.features:
           if feature.type == 'CDS':
               fid = feature_id(record, feature)
               print(fid)

is_gzip
~~~~~~~

.. code-block:: python

   is_gzip(gbkf)

Check if a file is gzip compressed.

**Parameters:**

- ``gbkf`` (str): Path to the file.

**Returns:**

- ``bool``: True if the file is gzip compressed, False otherwise.

**Example:**

.. code-block:: python

   from GenBankToLib.genbank import is_gzip
   
   if is_gzip('genome.gbk.gz'):
       print("File is compressed")

Data Structures
---------------

Bacteria Module
~~~~~~~~~~~~~~~

The ``bacteria`` module contains sets of Gram-positive and Gram-negative bacterial genera.

.. code-block:: python

   from GenBankToLib.bacteria import gram_positive, gram_negative
   
   # Check if a genus is Gram-positive
   if 'Escherichia' in gram_negative:
       print("Escherichia is Gram-negative")
   
   # Get all Gram-positive genera
   print(f"Gram-positive genera: {len(gram_positive)}")

Version Information
-------------------

.. code-block:: python

   from GenBankToLib import __version__
   
   print(f"genbank_to version: {__version__}")

Advanced Usage
--------------

Processing Large Files
~~~~~~~~~~~~~~~~~~~~~~

For large GenBank files, the generator-based functions are memory-efficient:

.. code-block:: python

   from GenBankToLib import genbank_to_faa
   
   # Process one protein at a time without loading entire file into memory
   protein_count = 0
   total_length = 0
   
   for seqid, protid, sequence in genbank_to_faa('large_genome.gbk'):
       protein_count += 1
       total_length += len(sequence)
   
   avg_length = total_length / protein_count
   print(f"Processed {protein_count} proteins, avg length: {avg_length:.1f}")

Custom Processing Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from GenBankToLib import genbank_to_faa, genbank_to_functions
   import json
   
   # Build a comprehensive protein database
   protein_db = {}
   
   # Get sequences
   for seqid, protid, sequence in genbank_to_faa('genome.gbk'):
       if protid not in protein_db:
           protein_db[protid] = {'sequence': str(sequence), 'length': len(sequence)}
   
   # Add functions
   for protid, function in genbank_to_functions('genome.gbk'):
       if protid in protein_db:
           protein_db[protid]['function'] = function
   
   # Save to JSON
   with open('protein_database.json', 'w') as f:
       json.dump(protein_db, f, indent=2)

Working with BioPython
~~~~~~~~~~~~~~~~~~~~~~

genbank_to functions integrate seamlessly with BioPython:

.. code-block:: python

   from GenBankToLib import genbank_seqio
   from Bio import SeqIO
   from Bio.SeqUtils import molecular_weight
   
   # Get molecular weights of all proteins
   parser, handle = genbank_seqio('genome.gbk')
   
   for record in parser:
       for feature in record.features:
           if feature.type == 'CDS' and 'translation' in feature.qualifiers:
               protein_seq = feature.qualifiers['translation'][0]
               mw = molecular_weight(protein_seq, seq_type='protein')
               print(f"{feature.qualifiers.get('protein_id', ['Unknown'])[0]}: {mw:.2f} Da")
   
   handle.close()

Error Handling
--------------

Example of proper error handling:

.. code-block:: python

   from GenBankToLib import genbank_to_faa
   import logging
   
   logging.basicConfig(level=logging.INFO)
   
   try:
       for seqid, protid, sequence in genbank_to_faa('genome.gbk'):
           # Process sequence
           pass
   except FileNotFoundError:
       logging.error("GenBank file not found")
   except Exception as e:
       logging.error(f"Error processing GenBank file: {e}")

Type Hints
----------

The library uses type hints for better IDE support:

.. code-block:: python

   from typing import Iterator, Tuple
   from Bio.Seq import Seq
   
   def process_proteins(gbkf: str) -> Iterator[Tuple[str, str, str]]:
       """
       Process proteins from a GenBank file.
       
       Args:
           gbkf: Path to GenBank file
           
       Yields:
           Tuple of (sequence_id, protein_id, sequence)
       """
       from GenBankToLib import genbank_to_faa
       return genbank_to_faa(gbkf)
