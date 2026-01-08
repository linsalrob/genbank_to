#!/usr/bin/env python
"""
Convert a genbank file to sequences
"""

import os
import sys
import gzip
import argparse
import json
from .genbank import genbank_to_faa, genbank_to_fna, genbank_to_orfs, genbank_to_ptt, genbank_to_functions
from .genbank import genbank_to_gff, genbank_to_phage_finder, genbank_seqio, genbank_to_amrfinder
from .genbank_to_json import genbank_to_json
from Bio import SeqIO
import logging
from .version import __version__

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


def run():
    """
    The entry point for an application
    """
    parser = argparse.ArgumentParser(description="Convert GenBank format files to...")
    parser.add_argument('-g', '--genbank', help='genbank file', required=True)
    parser.add_argument('-c', '--complex', help='complex identifier line', action='store_true')
    parser.add_argument('-a', '--aminoacids', help="output file for the amino acid sequences")
    parser.add_argument('-n', '--nucleotide', help='output file for nucleotide sequence')
    parser.add_argument('-p', '--ptt', help='output file for the ptt protein table')
    parser.add_argument('-o', '--orfs', help='output file for orfs')
    parser.add_argument('-f', '--functions', help='output file for two column table of [protein id, function]')
    parser.add_argument('-i', '--seqid', help='Only output these sequence ID(s) [multiple -i allowed]',
                        action='append')
    parser.add_argument('--gff3', help="Output gff3 format")
    parser.add_argument('--amr', help="Output NCBI AMRFinderPlus format (creates a gff file and an faa file)")
    parser.add_argument('--phage_finder', help='make a phage finder file')

    parser.add_argument('--bakta-json', help="Output BAKTA json format")
    bakta_group = parser.add_argument_group("BAKTA JSON metadata (requires --bakta-json)")
    # JSON-only arguments (Bakta metadata)
    bakta_group.add_argument('--bakta-version', default=None,
                        help='Bakta version string (default: NA) [requires --bakta-json]')
    bakta_group.add_argument('--db-version', default=None,
                        help='Database version string (default: NA) [requires --bakta-json]')

    bakta_group.add_argument('--genus', default=None,
                        help='Genus name (overrides GenBank annotation) [requires --bakta-json]')
    bakta_group.add_argument('--species', default=None,
                        help='Species name (overrides GenBank annotation) [requires --bakta-json]')
    bakta_group.add_argument('--strain', default=None,
                        help='Strain name (overrides GenBank annotation) [requires --bakta-json]')
    bakta_group.add_argument('--gram', default=None, choices=['+', '-'],
                        help='Gram stain (+ or -) [requires --bakta-json]')
    bakta_group.add_argument('--translation-table', type=int, default=None,
                        help='NCBI translation table number (default: 11) [requires --bakta-json]')


    parser.add_argument('--pseudo', help='include pseudo genes in the output. (This may cause biopython errors).',
                        action='store_true')
    parser.add_argument('--separate', action='store_true',
                        help='separate output into different files (with no other options just output gbk).')
    parser.add_argument('-z', '--zip', help='gzip compress the output. Experimental and may not work with everything!',
                        action='store_true')
    parser.add_argument('--log', help='Log file. Default = genbank_to.log', type=str, default='genbank_to.log')
    parser.add_argument('-d', '--debug', help='enable debugging', action='store_true')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    args = parser.parse_args()

    if not os.path.exists(args.genbank):
        sys.stderr.write(f"FATAL: {args.genbank} does not exist. Please check the file path and try again!")
        sys.exit(1)

    if args.seqid and not args.separate:
        sys.stderr.write("-i was provided, so requiring to separate files (--separate assumed)\n")
        args.separate = True

    json_only_fields = ['bakta_version', 'db_version', 'genus', 'species', 'strain', 'gram', 'translation_table']

    json_enabled = bool(args.bakta_json is not None)

    if not json_enabled:
        used = [f'--{name.replace("_", "-")}' for name in json_only_fields
                if getattr(args, name) is not None]
        if used:
            parser.error(f"{', '.join(used)} require --bakta-json")

    loglevel = logging.INFO
    if args.debug:
        loglevel = logging.DEBUG

    logging.basicConfig(filename=args.log, level=loglevel)
    logging.info(f"Starting genbank_to with {sys.argv}")

    did = False
    if args.nucleotide:
        if args.separate:
            lastid = None
            out = None
            for sid, seq in genbank_to_fna(args.genbank):
                if args.seqid and sid not in args.seqid:
                    continue
                if sid != lastid:
                    if out:
                        out.close()
                    out = open(f"{args.nucleotide}.{sid}.fna", 'w')
                    lastid = sid
                logging.info(f"Writing nucleotides for {sid} to {out.name}")
                out.write(f">{sid}\n{seq}\n")
            if out:
                out.close()
        else:
            with open(args.nucleotide, 'w') as out:
                logging.info(f"Writing all nucleotides to {args.nucleotide}")
                for sid, seq in genbank_to_fna(args.genbank):
                    out.write(f">{sid}\n{seq}\n")
        did = True

    if args.aminoacids:
        if args.separate:
            lastid = None
            out = None
            for seqid, sid, seq in genbank_to_faa(args.genbank, complexheader=args.complex, skip_pseudo=args.pseudo):
                if args.seqid and sid not in args.seqid:
                    continue
                if seqid != lastid:
                    if out:
                        out.close()
                    out = open(f"{args.aminoacids}.{seqid}.faa", 'w')
                    lastid = seqid
                logging.info(f"Writing amino acids for {sid} to {out.name}")
                out.write(f">{sid}\n{seq}\n")
            if out:
                out.close()
        else:
            with open(args.aminoacids, 'w') as out:
                logging.info(f"Writing all amino acids to {out.name}")
                for seqid, sid, seq in genbank_to_faa(args.genbank, complexheader=args.complex, skip_pseudo=args.pseudo):
                    out.write(f">{sid}\n{seq}\n")
        did = True

    if args.orfs:
        if args.separate:
            lastid = None
            out = None
            for seqid, sid, seq in genbank_to_orfs(args.genbank, complexheader=args.complex, skip_pseudo=args.pseudo):
                if args.seqid and sid not in args.seqid:
                    continue
                if seqid != lastid:
                    if out:
                        out.close()
                    out = open(f"{args.orfs}.{seqid}.orfs", 'w')
                    lastid = seqid
                logging.info(f"Writing orfs for {sid} to {out.name}")
                out.write(f">{sid}\n{seq}\n")
            if out:
                out.close()
        else:
            with open(args.orfs, 'w') as out:
                logging.info(f"Writing all orfs to {out.name}")
                for seqid, sid, seq in genbank_to_orfs(args.genbank, complexheader=args.complex, skip_pseudo=args.pseudo):
                    out.write(f">{sid}\n{seq}\n")
        did = True

    if args.ptt:
        r = genbank_to_ptt(args.genbank, False)
        logging.info(f"Writing ptt to {args.ptt}")
        with open(args.ptt, 'w') as out:
            for ln in r:
                out.write("\t".join(map(str, ln)))
                out.write("\n")
        did = True

    if args.functions:
        try:
            if args.zip:
                out = gzip.open(f"{args.functions}.gz", 'wt')
            else:
                out = open(args.functions, 'w')
            logging.info(f"Writing functions to {args.functions}")
            for sid, pid, prod in genbank_to_functions(args.genbank, True):
                out.write(f"{sid}\t{pid}\t{prod}\n")
            did = True
            out.close()
        except IOError as e:
            sys.stderr.write(f"There was an error writing to {args.functions}: {e}\n")
            sys.exit(1)

    if args.phage_finder:
        with open(args.phage_finder, 'w') as out:
            logging.info(f"Writing phage_finder to {args.phage_finder}")
            for tple in genbank_to_phage_finder(args.genbank):
                out.write("\t".join(map(str, tple)) + "\n")
        did = True

    if args.gff3:
        genbank_to_gff(args.genbank, args.gff3)
        did = True

    if args.bakta_json:
        bakta_json = genbank_to_json(args.genbank, args)
        with open(args.bakta_json, 'w') as out:
            logging.info(f"Writing bakta_json to {args.bakta_json}")
            json.dump(bakta_json, out, indent=4)
        did = True

    if args.amr:
        genbank_to_amrfinder(args.genbank, args.amr)
        did = True

    if not did and args.separate:
        for seq in genbank_seqio(args.genbank):
            if args.seqid and seq.id not in args.seqid:
                continue
            out = open(f"{seq.id}.gbk", 'w')
            logging.info(f"Writing {seq.id} to {out.name}")
            SeqIO.write(seq, out, 'genbank')
            out.close()
        did = True

    if not did:
        logging.warn("No option found\n")
        logging.warn("Please provide either a -n, -a, -o, -p, -f, --gff3 output file! (or all)")
        sys.exit(2)

    sys.exit(0)
