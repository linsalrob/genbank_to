Examples
========

This page provides comprehensive examples for various use cases of genbank_to.

Basic Examples
--------------

Convert GenBank to FASTA
~~~~~~~~~~~~~~~~~~~~~~~~~

The simplest use case - extract the genome sequence:

.. code-block:: bash

   genbank_to -g genome.gbk -n genome.fna

Extract All Sequence Types
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Get genome, ORFs, and proteins in one command:

.. code-block:: bash

   genbank_to -g genome.gbk \
       -n genome.fna \
       -o orfs.fna \
       -a proteins.faa

Working with Phage Genomes
---------------------------

Complete Phage Analysis
~~~~~~~~~~~~~~~~~~~~~~~

Extract all relevant information from a phage genome:

.. code-block:: bash

   # Download example phage genome (phiX174)
   wget https://raw.githubusercontent.com/linsalrob/genbank_to/main/test/NC_001417.gbk
   
   # Convert to multiple formats
   genbank_to -g NC_001417.gbk \
       -n NC_001417.fna \
       -a NC_001417.faa \
       -f NC_001417_functions.tsv \
       --gff3 NC_001417.gff3 \
       --phage_finder NC_001417.pf

Phage Annotation Comparison
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Compare annotations from different sources:

.. code-block:: bash

   # Extract functions from GenBank
   genbank_to -g phage1.gbk -f phage1_functions.tsv
   genbank_to -g phage2.gbk -f phage2_functions.tsv
   
   # Compare in Python
   import pandas as pd
   
   df1 = pd.read_csv('phage1_functions.tsv', sep='\t', 
                     names=['seqid', 'protid', 'function'])
   df2 = pd.read_csv('phage2_functions.tsv', sep='\t',
                     names=['seqid', 'protid', 'function'])
   
   # Merge and compare
   merged = df1.merge(df2, on='protid', suffixes=('_1', '_2'))
   print(merged[merged['function_1'] != merged['function_2']])

Bacterial Genome Analysis
--------------------------

Complete Bacterial Genome Conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   genbank_to -g ecoli.gbk \
       -n ecoli_genome.fna \
       -a ecoli_proteins.faa \
       -o ecoli_orfs.fna \
       -f ecoli_functions.tsv \
       --gff3 ecoli.gff3 \
       --bakta-json ecoli.json \
       --genus Escherichia \
       --species coli \
       --strain K-12 \
       --gram -

AMR Analysis Pipeline
~~~~~~~~~~~~~~~~~~~~~

Prepare files for antimicrobial resistance analysis:

.. code-block:: bash

   # Convert to AMRFinder format
   genbank_to -g bacteria.gbk --amr bacteria_amr
   
   # Run AMRFinderPlus (if installed)
   amrfinder \
       -n bacteria_amr.fna \
       -p bacteria_amr.faa \
       -g bacteria_amr.gff \
       -o amr_results.txt

Multi-Genome Processing
------------------------

Batch Processing
~~~~~~~~~~~~~~~~

Process multiple GenBank files:

.. code-block:: bash

   #!/bin/bash
   
   for gbk in genomes/*.gbk; do
       base=$(basename "$gbk" .gbk)
       echo "Processing $base..."
       
       genbank_to -g "$gbk" \
           -n "output/${base}.fna" \
           -a "output/${base}.faa" \
           --gff3 "output/${base}.gff3"
   done

Processing Multi-Record Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Handle GenBank files with multiple chromosomes or contigs:

.. code-block:: bash

   # Separate into individual files
   genbank_to -g multi_chromosome.gbk --separate -n chromosome
   
   # This creates:
   # chromosome.chr1.fna
   # chromosome.chr2.fna
   # etc.

Extract Specific Chromosomes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Extract only certain sequences from a multi-record file:

.. code-block:: bash

   genbank_to -g multi_chromosome.gbk \
       -i NC_000001 \
       -i NC_000002 \
       -n selected_chromosomes.fna

Python Library Usage
--------------------

Basic Library Import
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from GenBankToLib import (
       genbank_to_faa,
       genbank_to_fna,
       genbank_to_orfs,
       genbank_to_functions,
       genbank_to_gff
   )

Extract and Filter Proteins
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from GenBankToLib import genbank_to_faa
   
   # Extract proteins longer than 100 amino acids
   with open('large_proteins.faa', 'w') as out:
       for seqid, protid, sequence in genbank_to_faa('genome.gbk'):
           if len(sequence) > 100:
               out.write(f">{protid}\n{sequence}\n")

Build Custom Function Database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from GenBankToLib import genbank_to_functions
   import json
   
   # Build a dictionary of functions
   func_db = {}
   for protid, function in genbank_to_functions('genome.gbk'):
       func_db[protid] = function
   
   # Save as JSON
   with open('function_db.json', 'w') as f:
       json.dump(func_db, f, indent=2)

Calculate Genome Statistics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from GenBankToLib import genbank_to_fna, genbank_to_faa
   from Bio.SeqUtils import gc_fraction
   
   # Get genome statistics
   for seqid, sequence in genbank_to_fna('genome.gbk'):
       gc_content = gc_fraction(sequence) * 100
       length = len(sequence)
       print(f"{seqid}: {length} bp, {gc_content:.2f}% GC")
   
   # Count proteins
   protein_count = sum(1 for _ in genbank_to_faa('genome.gbk'))
   print(f"Proteins: {protein_count}")

Merge Multiple Genomes
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from GenBankToLib import genbank_to_faa
   
   genbank_files = ['genome1.gbk', 'genome2.gbk', 'genome3.gbk']
   
   with open('all_proteins.faa', 'w') as out:
       for gbk in genbank_files:
           for seqid, protid, sequence in genbank_to_faa(gbk):
               out.write(f">{protid}\n{sequence}\n")

Advanced Examples
-----------------

Complex Header Format
~~~~~~~~~~~~~~~~~~~~~

Generate detailed FASTA headers:

.. code-block:: bash

   genbank_to -g genome.gbk \
       -a proteins.faa \
       --complex

Output includes organism, location, product, and database IDs:

.. code-block:: text

   >NP_040703.1 [NC_001417] [Enterobacteria phage phiX174] [NC_001417_51_1905] 
   DNA replication protein [GeneID:1261050]

Include Pseudogenes
~~~~~~~~~~~~~~~~~~~

By default, pseudogenes are skipped. To include them:

.. code-block:: bash

   genbank_to -g genome.gbk -a proteins.faa --pseudo

Compressed Output
~~~~~~~~~~~~~~~~~

Generate compressed output (experimental):

.. code-block:: bash

   genbank_to -g genome.gbk \
       -f functions.tsv \
       --zip

Custom Logging
~~~~~~~~~~~~~~

Control logging output:

.. code-block:: bash

   # Custom log file with debug information
   genbank_to -g genome.gbk \
       -n genome.fna \
       --log my_conversion.log \
       --debug
   
   # View the log
   tail -f my_conversion.log

Integration Examples
--------------------

BLAST Database Creation
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Extract proteins
   genbank_to -g genome.gbk -a proteins.faa
   
   # Create BLAST database
   makeblastdb -in proteins.faa -dbtype prot -out genome_db

Genome Browser Setup
~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Generate GFF3 and FASTA
   genbank_to -g genome.gbk \
       -n genome.fna \
       --gff3 genome.gff3
   
   # Index for IGV
   samtools faidx genome.fna
   
   # Or sort GFF for JBrowse
   sort -k1,1 -k4,4n genome.gff3 > genome.sorted.gff3
   bgzip genome.sorted.gff3
   tabix -p gff genome.sorted.gff3.gz

Codon Usage Analysis
~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Extract ORFs
   genbank_to -g genome.gbk -o orfs.fna
   
   # Analyze with CodonW or similar tool
   codonw orfs.fna -all_indices

Phylogenetic Analysis
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Extract specific gene
   genbank_to -g genome1.gbk -a proteins1.faa
   genbank_to -g genome2.gbk -a proteins2.faa
   
   # Extract rpoB sequences (example)
   grep -A1 "RNA polymerase beta" proteins1.faa > rpoB.faa
   grep -A1 "RNA polymerase beta" proteins2.faa >> rpoB.faa
   
   # Align with MAFFT
   mafft rpoB.faa > rpoB_aligned.faa
   
   # Build tree with RAxML
   raxmlHPC -s rpoB_aligned.faa -n tree -m PROTGAMMAAUTO

Troubleshooting Examples
-------------------------

Check File Format
~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Verify it's a GenBank file
   head -1 genome.gbk
   # Should start with LOCUS
   
   # Check for multiple records
   grep -c "//" genome.gbk

Handle Gzipped Input
~~~~~~~~~~~~~~~~~~~~

genbank_to automatically handles gzipped files:

.. code-block:: bash

   genbank_to -g genome.gbk.gz -n genome.fna

Extract from Failed Conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If conversion fails partway through:

.. code-block:: bash

   # Enable debug logging
   genbank_to -g genome.gbk \
       -a proteins.faa \
       --debug \
       --log debug.log
   
   # Check the log for specific errors
   grep ERROR debug.log

Performance Optimization
------------------------

Large Genome Processing
~~~~~~~~~~~~~~~~~~~~~~~

For very large genomes or many files:

.. code-block:: bash

   # Process one output at a time to save memory
   genbank_to -g large_genome.gbk -n genome.fna
   genbank_to -g large_genome.gbk -a proteins.faa
   
   # Or use parallel processing for multiple files
   parallel -j 4 'genbank_to -g {} -n {.}.fna' ::: genomes/*.gbk

Streaming Processing
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from GenBankToLib import genbank_to_faa
   
   # Process proteins on-the-fly without loading all into memory
   for seqid, protid, sequence in genbank_to_faa('large_genome.gbk'):
       # Process each protein immediately
       result = analyze_protein(sequence)
       print(f"{protid}: {result}")
