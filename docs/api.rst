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
       genbank_seqio,
       genbank_to_json
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

genbank_to_json
~~~~~~~~~~~~~~~

.. code-block:: python

   genbank_to_json(genbank_path, genome_info)

Convert GenBank file to JSON format with comprehensive metadata.

**Parameters:**

- ``genbank_path`` (str): Path to the GenBank file.
- ``genome_info`` (dict): Dictionary containing optional genome information:
  
  - ``gram`` (str or None): Gram stain result ('+' for Gram-positive, '-' for Gram-negative, or None to infer)
  - ``translation_table`` (int or None): NCBI translation table number (default: 11)

**Returns:**

- ``dict``: Complete JSON data structure containing:
  
  - ``genome``: Genome metadata (genus, species, strain, gram, translation_table, complete)
  - ``stats``: Genome statistics (no_sequences, size, gc, n_ratio, n50, coding_ratio)
  - ``features``: List of feature dictionaries (CDS, tRNA, rRNA, tmRNA, ncRNA)
  - ``sequences``: List of sequence objects with metadata
  - ``version``: Version information

**Example:**

.. code-block:: python

   from GenBankToLib import genbank_to_json
   import json
   
   # Basic usage with defaults
   genome_info = {
       'gram': None,
       'translation_table': 11
   }
   json_data = genbank_to_json('genome.gbk', genome_info)
   
   # Save to file
   with open('genome.json', 'w') as f:
       json.dump(json_data, f, indent=2)
   
   # Specify Gram stain and translation table
   genome_info = {
       'gram': '-',
       'translation_table': 11
   }
   json_data = genbank_to_json('ecoli.gbk', genome_info)
   
   # Access the data
   print(f"Genome: {json_data['genome']['genus']} {json_data['genome']['species']}")
   print(f"Size: {json_data['stats']['size']} bp")
   print(f"GC content: {json_data['stats']['gc']:.2%}")
   print(f"Features: {len(json_data['features'])}")

**JSON Output Structure:**

The output JSON contains the following structure:

.. code-block:: json

   {
       "genome": {
           "genus": "Escherichia",
           "species": "coli",
           "strain": "K-12",
           "complete": true,
           "gram": "-",
           "translation_table": 11
       },
       "stats": {
           "no_sequences": 1,
           "size": 4641652,
           "gc": 0.5079,
           "n_ratio": 0.0,
           "n50": 4641652,
           "coding_ratio": 0.8756
       },
       "features": [
           {
               "type": "CDS",
               "contig": "NC_000913",
               "start": 190,
               "stop": 255,
               "strand": 1,
               "gene": "thrL",
               "locus": "b0001",
               "product": "thr operon leader peptide",
               "nt": "ATGAAACGC...",
               "aa": "MKRISTT...",
               "aa_hexdigest": "abc123...",
               "frame": 0,
               "start_type": "ATG",
               "rbs_motif": "AGGAG"
           }
       ],
       "sequences": [
           {
               "id": "NC_000913",
               "description": "Escherichia coli str. K-12 substr. MG1655, complete genome",
               "sequence": "AGCTTTTCATTC...",
               "length": 4641652,
               "complete": true,
               "type": "chromosome",
               "topology": "circular"
           }
       ],
       "version": {
           "genbank_to": "1.2.3"
       }
   }

**Notes:**

- The JSON format is compatible with Bakta annotation output
- Coordinates are 1-based inclusive (start, stop)
- Frame is 0-based (0, 1, 2) where GenBank codon_start values 1, 2, 3 are converted to frame values 0, 1, 2 respectively
- GC content is calculated as (G+C)/(A+C+G+T), ignoring ambiguous bases
- MD5 hexdigest is provided for amino acid sequences in the ``aa_hexdigest`` field
- Gram stain can be inferred from genus if not provided

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

   from GenBankToLib import genbank_to_json
   import json
   
   # Build a comprehensive genome database with genbank_to_json
   genome_info = {
       'gram': '-',
       'translation_table': 11
   }
   
   # Convert to JSON format
   genome_data = genbank_to_json('genome.gbk', genome_info)
   
   # Process the data
   print(f"Genome: {genome_data['genome']['genus']} {genome_data['genome']['species']}")
   print(f"Total size: {genome_data['stats']['size']} bp")
   print(f"GC content: {genome_data['stats']['gc']:.2%}")
   print(f"Coding ratio: {genome_data['stats']['coding_ratio']:.2%}")
   
   # Extract specific features
   cds_features = [f for f in genome_data['features'] if f['type'] == 'CDS']
   print(f"CDS features: {len(cds_features)}")
   
   # Save to JSON file
   with open('genome_analysis.json', 'w') as f:
       json.dump(genome_data, f, indent=2)

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
