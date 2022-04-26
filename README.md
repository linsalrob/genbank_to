# genbank_to

[![Edwards Lab](https://img.shields.io/badge/Bioinformatics-EdwardsLab-03A9F4)](https://edwards.flinders.edu.au/)
[![DOI](https://www.zenodo.org/badge/481464683.svg)](https://www.zenodo.org/badge/latestdoi/481464683)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![GitHub language count](https://img.shields.io/github/languages/count/linsalrob/genbank_to)
[![PyPi](https://img.shields.io/pypi/pyversions/genbank_to.svg?style=flat-square&label=PyPi%20Versions)](https://pypi.org/project/genbank_to/)

A straightforward application to convert NCBI GenBank format files to a swath of other formats. Hopefully we have the 
format you need, but if not either post [an issue](https://github.com/linsalrob/genbank_to/issues) using our template,
or if you have already got it working, post [a PR](https://github.com/linsalrob/genbank_to/pulls) so we can add it and
add you to the project.

You might also be interested [deprekate's](https://github.com/deprekate/) package called [genbank](https://github.com/deprekate/genbank) which includes
several of the features here, and you can `import genbank` into your Python projects.

# What it does

Read an NCBI GenBank format file (like our [test data](test/NC_001417.gbk)) and convert it to one of many
different formats.

# Input formats

At the moment we only support NCBI GenBank format. If you want us to read other common formats, 
[let us know](https://github.com/linsalrob/genbank_to/issues) and we'll add them.

# Output formats

Here are the output formats you can request. You can request as many of these at once as you like!

These outputs are assuming you provide a (for example) genome file that contains ORFs, Proteins, and Genomes.

## Nucleotide output

 - `-n` or `--nucleotide` outputs the whole DNA sequence (e.g. the genome)
 - `-o` or `--orfs` outputs the DNA sequence of the open reading frames

## Protein output

 - `-a` or `--aminoacids` outputs the protein sequence for each of the open reading frames

## Complex formats

 - `-p` or `--ptt` NCBI ptt protein table. This is a somewhat deprecated NCBI format from their genomes downloads
 - `-f` or `--functions` outputs tab separated data of `protein ID` and `protein function` (also called the `product`)
 - `--gff3` outputs GFF3 format
 - `--amr`  outputs a GFF file, an amino acid fasta file, and a nucleotide fasta file as required by [AMR Finder Plus](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#examples). Note that this format checks for validity that often crashes AMRFinderPlus
 - `--phage_finder` outputs a unique format required by [phage_finder](http://phage-finder.sourceforge.net/)

## Output options

 - `--pseudo` normally we skip pseudogenes (e.g. in creating amino acid fasta files). This will try and include pseudogenes, but often biopython complains and ignores them!
 - `-i` or `--seqid` only output this sequence, or these sequences if you specify more than one `-i`/`--seqid`
 - `-z` or `--zip` compress some of the outputs
 - `--log` write logs to a different file

## Separate multi-GenBank files

If your GenBank files contains multiple sequence records (separated with `//`), you can provide the `--separate` flag. 
This will write each entry into its own file. This is compatible with `-n`/`--nucleotide`, `-o`/`--orfs`, and
`-a`/`--aminoacids`. However, if you provide the `--separate` flag on its own, it will write each entry in your 
multi-GenBank file to its own GenBank file.

## Examples

All of these examples use our [test data](test/NC_001417.gbk)

1. Extract a `fasta` of the genome:

```bash
genbank_to -g test/NC_001417.gbk -n test/NC_001417.fna
```

2. Extract the DNA sequences of the ORFs to a single file

```bash
genbank_to -g test/NC_001417.gbk -o test/NC_001417.orfs
```

3. Extract the protein (amino acid) sequences of the ORFs to a file

```bash
genbank_to -g test/NC_001417.gbk -a test/NC_001417.faa
```

4. Do all of these at once

```bash
genbank_to -g test/NC_001417.gbk -n test/NC_001417.fna -o test/NC_001417.orfs -a test/NC_001417.faa
```

# Installation

You can install `genbank_to` in three different ways:

1. Using conda

This is the easiest and recommended method.

```bash
mamba create -n genbank_to genbank_to
conda activate genbank_to
genbank_to --help
```

2. Using pip

I recommend putting this into a virtual environment:

```bash
virtualenv venv
source venv/bin/activate
pip install genbank_to
genbank_to --help
```

3. Directly from this repository

(Not really recommended as things might break)

```bash
git clone https://github.com/linsalrob/genbank_to.git
cd genbank_to
virtualenv venv
source venv/bin/activate
python setup.py install
genbank_to --help
```


