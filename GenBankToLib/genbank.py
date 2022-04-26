"""
Read a genbank file and do things with it!
"""
from typing import Any

import sys
import re
import binascii
import gzip
from Bio import SeqIO
from BCBio import GFF  # bcbio-gff package
import pandas as pd
from io import StringIO
import logging

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


def is_gzip(gbkf):
    """
    Is the file compressed?
    :param gbkf:
    :return: true if compressed else false
    """

    t=""
    with open(gbkf, 'rb') as i:
        t = i.read(2)

    return binascii.hexlify(t) == b'1f8b'

def cds_details(seq, feat, complexheader=False, skip_pseudo=True):
    """
    Extract the feature details, and return them.

    Returns None if not a valid feature, so check return
    """
    if feat.type != 'CDS':
        return None

    cid = feature_id(seq, feat)

    if skip_pseudo and 'pseudo' in feat.qualifiers:
        # sys.stderr.write(f"Skipped pseudogene in {cid}\n")
        return None

    if complexheader:
        (start, stop, strand) = (feat.location.start.position, feat.location.end.position, feat.strand)

        loc = f"{start}_{stop}"
        if strand < 0:
            loc = f"{stop}_{start}"

        cid += f' [{seq.id}] '
        if 'organism' in seq.annotations:
            cid += f' [{seq.annotations["organism"]}]'
        cid += f' [{seq.id}_{loc}]'
        if 'product' in feat.qualifiers:
            cid += f' {feat.qualifiers["product"][0]}'
        else:
            cid += f' [hypothetical protein]'

    return cid

def genbank_seqio(gbkf):
    """
    Get the parser stream
    :param gbkf: genbank file
    :return: the genbank parser and the file handle (to close it)
    """

    if is_gzip(gbkf):
        handle = gzip.open(gbkf, 'rt')
    else:
        handle = open(gbkf, 'r')
    return SeqIO.parse(handle, "genbank"), handle


def feature_id(seq, feat):
    """
    Choose the appropriate id for the feature
    :param seq: the sequence
    :param feat: the feature
    :return: the id
    """

    if 'protein_id' in feat.qualifiers:
        return '|'.join(feat.qualifiers['protein_id'])
    elif 'locus_tag' in feat.qualifiers:
        return "|".join(feat.qualifiers['locus_tag'])
    elif 'db_xref' in feat.qualifiers:
        return '|'.join(feat.qualifiers['db_xref'])
    else:
        return seq.id + "." + str(feat.location)


def feat_to_text(feat, qual):
    if qual in feat.qualifiers:
        return " ".join(feat.qualifiers[qual])
    return "-"


def genbank_to_fna(gbkf, include_definition=False):
    """
    Parse a genbank file
    :param gbkf: genbank file
    :param include_definition: include the genbank definition line with the accession
    :return: a dict of the sequences
    """

    seqs, handle = genbank_seqio(gbkf)
    for seq in seqs:
        myid = seq.id
        if include_definition:
            myid += " " + seq.description
        yield myid, seq.seq
    handle.close()

def genbank_to_faa(gbkf, complexheader=False, skip_pseudo=True):
    """
    Parse a genbank file
    :param gbkf: the genbank file
    :param complexheader: more detail in the header
    :return: yield the protein id and sequence
    """

    seqs, handle = genbank_seqio(gbkf)
    for seq in seqs:
        for feat in seq.features:
            cid = cds_details(seq, feat, complexheader, skip_pseudo)
            if not cid:
                continue
            
            if 'translation' in feat.qualifiers:
                yield seq.id, cid, feat.qualifiers['translation'][0]
            elif len(feat.extract(seq).seq) % 3 != 0:
                logging.debug(f"WARNING: Length of {seq.id} --> {cid} is not a multiple of 3 (it is {len(feat.extract(seq).seq)})")
            else:
                yield seq.id, cid, str(feat.extract(seq).translate().seq)
    handle.close()

def genbank_to_functions(gbkf, seqid=False, skip_pseudo=True):
    """
    Parse a genbank file
    :param gbkf: the genbank file
    :param seqid: include the sequence ID in the yield
    :return: yield a tple of [(seqid), protein id, function]
    """
    seqs, handle = genbank_seqio(gbkf)
    for seq in seqs:
        for feat in seq.features:
            cid = cds_details(seq, feat, False, skip_pseudo)
            if not cid:
                continue

            prod = "Hypothetical protein"
            if "product" in feat.qualifiers:
                prod = "|".join(feat.qualifiers['product'])

            if seqid:
                yield seq.id, cid, prod
            else:
                yield cid, prod
    handle.close()

def genbank_to_orfs(gbkf, complexheader=False, skip_pseudo=True):
    """
    Parse a genbank file
    :param gbkf:
    :param complexheader:
    :return: a dict of the sequences
    """

    seqs, handle = genbank_seqio(gbkf)
    for seq in seqs:
        for feat in seq.features:
            cid = cds_details(seq, feat, complexheader, skip_pseudo)
            if not cid:
                continue

            yield seq.id, cid, str(feat.extract(seq).seq)
    handle.close()

def genbank_to_ptt(gbkf, printout=False):
    """
    Convert the genbank file to a table with the same columns of the ptt file
    :param gbkf: the genbank input file
    :param printout: print the table
    :return: the table
    """

    res = []

    gire = re.compile('GI:(\\d+)')
    cogre = re.compile('(COG\\S+)')
    seqs, handle = genbank_seqio(gbkf)
    for seq in seqs:
        for feat in seq.features:
            if feat.type != "CDS":
                continue

            gi = "-"
            if gire.match(feat_to_text(feat, 'db_xref')):
                gi = gire.match(feat_to_text(feat, 'db_xref'))[1]

            cog = "-"
            if cogre.match(feat_to_text(feat, 'product')):
                cog = cogre.match(feat_to_text(feat, 'product'))[1]

            gene = feat_to_text(feat, 'gene')
            if gene == "-":
                gene = str(feat.location)

            cid = feature_id(seq, feat)

            thisres = [
                f"{feat.location.start}..{feat.location.end}",
                "+" if feat.strand >= 0 else "-",
                (len(feat.location) / 3) - 1,
                gi,
                gene,
                cid,
                cog,
                feat_to_text(feat, 'product')
            ]

            if printout:
                print("\t".join(map(str, thisres)))

            res.append(thisres)

    handle.close()

    return res


def genbank_to_phage_finder(gbkf):
    """
    This is a very specific format used by phage_finder (https://phage-finder.sourceforge.net/documentation.htm)
    - contig id, including >
    - length of the contig
    - gene ID
    - start
    - end
    - function

    :param gbkf: the genbank file
    :return: yields a tple of this data
    """

    seqs, handle = genbank_seqio(gbkf)
    for seq in seqs:
        for feat in seq.features:
            if feat.type != 'CDS':
                continue
            cid = feature_id(seq, feat)
            fn = "Hypothetical protein"
            if 'product' in feat.qualifiers:
                fn = feat_to_text(feat, 'product')
            yield [seq.id, len(seq.seq), cid, feat.location.start, feat.location.end, fn]
    handle.close()

def genbank_to_pandas(gbkf, mincontiglen, ignorepartials=True, convert_selenocysteine=False):
    """
    This is a bit of a specific format used by phage_boost. its a simple dataframe with a couple of
    additional columns:
        ['contig',
         'id',
         'start',
         'stop',
         'direction',
         'partial',
         'DNAseq',
         'AAseq',
         'header']
    :param ignorepartials: Ignore any gene call with a frameshift (ie. a stop codon in the middle of the sequence)
    :param convert_selenocysteine: PhageBoost crashes with a selenocysteine protein because it is not included in
    Biopython
    :param gbkf: Genbank file to parse
    :param mincontiglen: minimum contig length to include
    :return: a pandas data frame
    """

    c = 0
    genes = []
    seqs, handle = genbank_seqio(gbkf)
    for seq in seqs:
        if len(seq) < mincontiglen:
            sys.stderr.write(f"Skipped {seq.id} because it's length ({len(seq)}) is less than "
                             f"the minimum contig length ({mincontiglen})\n")
            continue
        for feat in seq.features:
            if feat.type != 'CDS':
                continue

            tid = seq.id + "_" + str(c)
            partial = 0
            # I don't think this is exactly right
            if 'truncated' in feat.qualifiers:
                partial = 1

            dnaseq = str(feat.extract(seq).seq)
            if len(dnaseq) == 0:
                sys.stderr.write(f"The DNA sequence for {feature_id(seq, feat)} was zero, so skipped\n")
                continue

            # we just do a de novo translation rather than relying on the translation provided
            # in the genbank file that is often wrong
            trans = str(feat.extract(seq).translate().seq)

            while trans.endswith('*'):
                trans = trans[:-1]

            # partial amino acid codes we should ignore
            paa = {'B', 'Z', 'J', 'X', '*'}

            keeporf = True

            if ignorepartials:
                for aa in paa:
                    if aa in trans:
                        sys.stderr.write(f"There is a {aa} in  {feature_id(seq, feat)} so skipped.\n")
                        keeporf = False

            if not keeporf:
                continue

            if len(trans) == 0:
                sys.stderr.write(f"The translation for {feature_id(seq, feat)} was zero, so skipped\n")
                continue

            if convert_selenocysteine:
                trans = trans.replace('U', 'C')
            row = [seq.id, c, feat.location.start.position, feat.location.end.position, feat.strand,
                   partial, dnaseq, trans, tid]
            c += 1

            genes.append(row)

    handle.close()

    genecalls = pd.DataFrame(genes, columns=['contig', 'id', 'start', 'stop', 'direction', 'partial', 'DNAseq', 'AAseq',
                                             'header'])

    return genecalls


def genbank_to_gff(gbkf, out_gff):
    """
    Convert the genbank file to a gff3 format file
    :param gbkf: the genbank input file
    :param out_gff: the gff3 format output file
    :param verbose: more output
    :return: nothing
    """

    logging.info(f"Writing gff3 to {out_gff}")
    with open(out_gff, 'w') as outf:
        seqs, handle = genbank_seqio(gbkf)
        GFF.write(seqs, outf, True)
        handle.close()


def genbank_to_amrfinder(gbkf, amrout):
    """
    Convert the genbank file to amr finder format. This is a bastardized GFF3 format that does
    not include the ##FASTA header, and also requires all the proteins to be present
    and have their protein IDs as a "Name" field (caps seem important) in the GFF
    """

    logging.info(f"AMR: Converting {gbkf} to GFF3")
    gffstr = StringIO()
    seqs, handle = genbank_seqio(gbkf)
    GFF.write(seqs, gffstr, True)
    handle.close()
    gffstr.seek(0)

    logging.info(f"AMR: Converting {gbkf} to amino acids in {amrout}.faa")
    seqs = set()
    with open(f"{amrout}.faa", 'w') as out:
        for seqid, sid, seq in genbank_to_faa(gbkf, False, skip_pseudo=True):
            seqs.add(sid)
            out.write(f">{sid}\n{seq}\n")

    logging.info(f"AMR: Converting {gbkf} to nucleotides in {amrout}.fna")
    with open(f"{amrout}.fna", 'w') as out:
        for sid, seq in genbank_to_fna(gbkf):
            out.write(f">{sid}\n{seq}\n")

    logging.info(f"Writing the GFF to {amrout}.gff")
    search = re.compile(r'protein_id=([\w\.]+)')
    with open(f"{amrout}.gff", 'w') as out:
        for r in gffstr:
            if r.startswith("##FASTA"):
                break
            p: List = r.split("\t")
            if len(p) > 1 and p[2] == 'CDS':
                m = search.search(p[8])
                if m and m.groups()[0] in seqs:
                    tmp = p[8].rstrip()
                    tmp += f";Name={m.groups()[0]}\n"
                    p[8] = tmp
                    out.write("\t".join(map(str, p)))
                elif m:
                    logging.debug(f"Could not find {m.groups()[0]} in {seqs}")
                elif 'pseudo=' in p[8]:
                    pass
                else:
                    logging.debug(f"Could not parse protein_id from {p[8]}")
            else:
                out.write(r)
