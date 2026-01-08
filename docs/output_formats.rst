Output Formats
==============

This page describes all the output formats supported by genbank_to in detail.

FASTA Formats
-------------

Nucleotide FASTA (Genome)
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Option**: ``-n``, ``--nucleotide``

**Description**: Outputs the complete nucleotide sequence(s) from the GenBank file.

**Format**: Standard FASTA format with sequence ID from the GenBank LOCUS or ACCESSION field.

**Example Output**:

.. code-block:: text

   >NC_001417.1
   GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT
   GATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGACTGCTGGCGGAAAATGAGAAA
   ...

**Use Cases**:

- Reference genome for mapping
- Input for genome assembly comparison
- Sequence extraction for primers/probes

ORF Nucleotide FASTA
~~~~~~~~~~~~~~~~~~~~

**Option**: ``-o``, ``--orfs``

**Description**: Outputs the DNA sequences of all coding sequences (CDS features).

**Format**: FASTA format with protein IDs as headers (or locus tags if protein IDs are not available).

**Example Output**:

.. code-block:: text

   >NP_040703.1
   ATGGTTAGCAAAATCGAACGTGCAAAGATTGATGATATTAATATTTTTATTGAAAATCACCAGAAAGATA
   TAGACTATCTTTGGCAACGTATACCGATGAAATCATTAAAGACTTAAAAGTTGAGCGCTTTGATACGAGT
   ...

**Use Cases**:

- Gene cloning
- Codon usage analysis
- Primer design
- Gene synthesis

Protein FASTA
~~~~~~~~~~~~~

**Option**: ``-a``, ``--aminoacids``

**Description**: Outputs amino acid sequences for all coding sequences.

**Format**: FASTA format with protein IDs. If ``--complex`` is specified, includes detailed annotation.

**Example Output (Simple)**:

.. code-block:: text

   >NP_040703.1
   MVSKIERCKILMIINFLIEIHQKDIDYLWQRIPEIIIKDLKVERFDDTVKVGGYKKGGLVQPGGSLRLYE
   VDEKGHFPENVVYDGDTVVADDTLYLVAVLDERKMKGINTRELLESYFDRRGFRLPVGHIDNKPGFNVK
   *

**Example Output (Complex with** ``--complex`` **)**:

.. code-block:: text

   >NP_040703.1 [NC_001417] [Enterobacteria phage phiX174] [NC_001417_51_1905] 
   gpA DNA replication protein [GeneID:1261050]

**Use Cases**:

- Protein homology searches (BLAST)
- Phylogenetic analysis
- Functional annotation
- Structure prediction

Structured Formats
------------------

GFF3 (Generic Feature Format)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Option**: ``--gff3``

**Description**: Outputs annotations in GFF3 format, a standardized format for genomic features.

**Format**: Tab-delimited with 9 columns: seqid, source, type, start, end, score, strand, phase, attributes.

**Example Output**:

.. code-block:: text

   ##gff-version 3
   ##sequence-region NC_001417 1 5386
   NC_001417	GenBank	region	1	5386	.	+	.	ID=NC_001417:1..5386;...
   NC_001417	GenBank	gene	51	1905	.	+	.	ID=gene-phiX174p01;...
   NC_001417	GenBank	CDS	51	1905	.	+	0	ID=cds-NP_040703.1;...

**Use Cases**:

- Genome browsers (IGV, JBrowse)
- Comparative genomics tools
- Annotation pipelines
- Data interchange

NCBI PTT Format
~~~~~~~~~~~~~~~

**Option**: ``-p``, ``--ptt``

**Description**: Protein table format formerly used by NCBI for genome downloads.

**Format**: Tab-delimited table with columns: Location, Strand, Length, PID, Gene, Synonym, COG, Product.

**Example Output**:

.. code-block:: text

   51..1905	+	617	-	gpA	NP_040703.1	-	DNA replication protein
   1906..2079	+	57	-	gpB	NP_040704.1	-	capsid morphogenesis protein
   2092..2529	+	145	-	gpD	NP_040705.1	-	capsid morphogenesis protein

**Use Cases**:

- Legacy pipeline compatibility
- Quick protein overview
- Tab-delimited data processing

Function Table
~~~~~~~~~~~~~~

**Option**: ``-f``, ``--functions``

**Description**: Simple two-column table mapping protein IDs to their functional annotations.

**Format**: Tab-separated values with columns: Protein ID, Function.

**Example Output**:

.. code-block:: text

   NC_001417	NP_040703.1	DNA replication protein
   NC_001417	NP_040704.1	capsid morphogenesis protein
   NC_001417	NP_040705.1	capsid morphogenesis protein
   NC_001417	NP_040706.1	DNA maturase protein B

**Use Cases**:

- Functional enrichment analysis
- Database imports
- Quick annotation lookup
- Spreadsheet analysis

Specialized Formats
-------------------

Bakta JSON
~~~~~~~~~~

**Option**: ``--bakta-json``

**Description**: JSON format compatible with Bakta genome annotation output. Includes comprehensive metadata and feature annotations.

**Additional Options**:

- ``--bakta-version``: Version string
- ``--db-version``: Database version
- ``--genus``, ``--species``, ``--strain``: Organism information
- ``--gram``: Gram stain (+/-)
- ``--translation-table``: Genetic code

**Example Output**:

.. code-block:: json

   {
       "version": "1.0",
       "genome": {
           "genus": "Enterobacteria",
           "species": "phage phiX174",
           "strain": "Sangier",
           "gram": "-"
       },
       "sequences": [
           {
               "id": "NC_001417",
               "length": 5386,
               "gc": 0.447,
               "features": [
                   {
                       "type": "cds",
                       "start": 51,
                       "stop": 1905,
                       "strand": "+",
                       "product": "DNA replication protein"
                   }
               ]
           }
       ]
   }

**Use Cases**:

- Bakta pipeline integration
- Structured data analysis
- Web applications
- Database storage

AMRFinderPlus Format
~~~~~~~~~~~~~~~~~~~~

**Option**: ``--amr``

**Description**: Creates three files required by NCBI's AMRFinderPlus tool for antimicrobial resistance gene annotation.

**Output Files**:

1. ``BASENAME.gff`` - Modified GFF3 with Name attributes
2. ``BASENAME.faa`` - Protein sequences
3. ``BASENAME.fna`` - Nucleotide sequences

**Special Features**:

- Validates format for AMRFinderPlus compatibility
- Adds required Name fields to GFF
- Excludes pseudogenes

**Use Cases**:

- AMR gene detection
- Resistance profile analysis
- Public health surveillance
- Clinical microbiology

Phage Finder Format
~~~~~~~~~~~~~~~~~~~

**Option**: ``--phage_finder``

**Description**: Tab-delimited format required by the phage_finder tool for prophage identification.

**Format**: Tab-separated with columns: Contig ID, Contig Length, Gene ID, Start, End, Function.

**Example Output**:

.. code-block:: text

   NC_001417	5386	NP_040703.1	51	1905	DNA replication protein
   NC_001417	5386	NP_040704.1	1906	2079	capsid morphogenesis protein
   NC_001417	5386	NP_040705.1	2092	2529	capsid morphogenesis protein

**Use Cases**:

- Prophage detection in bacterial genomes
- Phage-host interaction studies
- Comparative phage genomics

Output Modifiers
----------------

Separate Files
~~~~~~~~~~~~~~

**Option**: ``--separate``

**Description**: When working with multi-record GenBank files, creates separate output files for each sequence record.

**Naming Convention**: ``BASENAME.SEQID.EXTENSION``

**Example**:

.. code-block:: bash

   genbank_to -g multi.gbk --separate -n output
   # Creates: output.NC_001417.fna, output.NC_001418.fna, etc.

Sequence ID Filtering
~~~~~~~~~~~~~~~~~~~~~

**Option**: ``-i``, ``--seqid``

**Description**: Filters output to include only specified sequence IDs. Can be used multiple times.

**Example**:

.. code-block:: bash

   genbank_to -g multi.gbk -i NC_001417 -i NC_001418 -n output.fna

Complex Headers
~~~~~~~~~~~~~~~

**Option**: ``--complex``

**Description**: Adds detailed information to FASTA headers including organism, location, product, and database cross-references.

**Example**:

.. code-block:: text

   >NP_040703.1 [NC_001417] [Enterobacteria phage phiX174] [NC_001417_51_1905] 
   DNA replication protein [GeneID:1261050]

Compression
~~~~~~~~~~~

**Option**: ``-z``, ``--zip``

**Description**: Compresses output files using gzip. Experimental feature.

**Example**:

.. code-block:: bash

   genbank_to -g genome.gbk -f functions.tsv --zip
   # Creates: functions.tsv.gz

Format Comparison
-----------------

.. list-table:: Output Format Comparison
   :header-rows: 1
   :widths: 20 20 20 40

   * - Format
     - Type
     - Primary Use
     - Tools
   * - Nucleotide FASTA
     - Sequence
     - Genome storage
     - BWA, Bowtie, BLAST
   * - Protein FASTA
     - Sequence
     - Protein analysis
     - BLAST, DIAMOND, InterProScan
   * - ORF FASTA
     - Sequence
     - Gene analysis
     - Gene synthesis, primers
   * - GFF3
     - Annotation
     - Feature storage
     - IGV, JBrowse, bedtools
   * - PTT
     - Table
     - Legacy compatibility
     - Custom scripts
   * - Functions
     - Table
     - Annotation lookup
     - Spreadsheets, R/Python
   * - Bakta JSON
     - Structured
     - Data interchange
     - Bakta, web apps
   * - AMRFinder
     - Specialized
     - AMR detection
     - AMRFinderPlus
   * - Phage Finder
     - Specialized
     - Prophage detection
     - phage_finder
