#!/usr/bin/env python

##### findSpanin.pl --> findSpanin.py
######### Incooperated from the findSpanin.pl script, but better and more snakey.

import argparse
from cpt import OrfFinder
from Bio import SeqIO
from Bio import Seq
import re
from statistics import median
from spaninFuncs import *
import os

# if __name__ == '__main__':
# pass
###############################################################################

if __name__ == "__main__":

    # Common parameters for both ISP / OSP portion of script

    parser = argparse.ArgumentParser(
        description="Get putative protein candidates for spanins"
    )
    parser.add_argument(
        "fasta_file", type=argparse.FileType("r"), help="Fasta file"
    )  # the "input" argument

    parser.add_argument(
        "-f",
        "--format",
        dest="seq_format",
        default="fasta",
        help="Sequence format (e.g. fasta, fastq, sff)",
    )  # optional formats for input, currently just going to do ntFASTA

    parser.add_argument(
        "--strand",
        dest="strand",
        choices=("both", "forward", "reverse"),
        default="both",
        help="select strand",
    )  # Selection of +, -, or both strands

    parser.add_argument(
        "--table", dest="table", default=11, help="NCBI Translation table", type=int
    )  # Uses "default" NCBI codon table. This should always (afaik) be what we want...

    parser.add_argument(
        "-t",
        "--ftype",
        dest="ftype",
        choices=("CDS", "ORF"),
        default="ORF",
        help="Find ORF or CDSs",
    )  # "functional type(?)" --> Finds ORF or CDS, for this we want just the ORF

    parser.add_argument(
        "-e",
        "--ends",
        dest="ends",
        choices=("open", "closed"),
        default="closed",
        help="Open or closed. Closed ensures start/stop codons are present",
    )  # includes the start and stop codon

    parser.add_argument(
        "-m",
        "--mode",
        dest="mode",
        choices=("all", "top", "one"),
        default="all",  # I think we want this to JUST be all...nearly always
        help="Output all ORFs/CDSs from sequence, all ORFs/CDSs with max length, or first with maximum length",
    )

    parser.add_argument(
        "--switch",
        dest="switch",
        default="all",
        help="switch between ALL putative osps, or a range. If not all, insert a range of two integers separated by a colon (:). Eg: 1234:4321",
    )
    # isp parameters
    parser.add_argument(
        "--isp_min_len",
        dest="isp_min_len",
        default=60,
        help="Minimum ORF length, measured in codons",
        type=int,
    )
    parser.add_argument(
        "--isp_on",
        dest="out_isp_nuc",
        type=argparse.FileType("w"),
        default="_out_isp.fna",
        help="Output nucleotide sequences, FASTA",
    )
    parser.add_argument(
        "--isp_op",
        dest="out_isp_prot",
        type=argparse.FileType("w"),
        default="_out_isp.fa",
        help="Output protein sequences, FASTA",
    )
    parser.add_argument(
        "--isp_ob",
        dest="out_isp_bed",
        type=argparse.FileType("w"),
        default="_out_isp.bed",
        help="Output BED file",
    )
    parser.add_argument(
        "--isp_og",
        dest="out_isp_gff3",
        type=argparse.FileType("w"),
        default="_out_isp.gff3",
        help="Output GFF3 file",
    )
    parser.add_argument(
        "--isp_min_dist",
        dest="isp_min_dist",
        default=10,
        help="Minimal distance to first AA of TMD, measured in AA",
        type=int,
    )
    parser.add_argument(
        "--isp_max_dist",
        dest="isp_max_dist",
        default=30,
        help="Maximum distance to first AA of TMD, measured in AA",
        type=int,
    )
    parser.add_argument(
        "--putative_isp",
        dest="putative_isp_fa",
        type=argparse.FileType("w"),
        default="_putative_isp.fa",
        help="Output of putative FASTA file",
    )
    parser.add_argument(
        "--min_tmd_size",
        dest="min_tmd_size",
        default=10,
        help="Minimal size of the TMD domain",
        type=int,
    )
    parser.add_argument(
        "--max_tmd_size",
        dest="max_tmd_size",
        default=20,
        help="Maximum size of the TMD domain",
        type=int,
    )
    parser.add_argument(
        "--summary_isp_txt",
        dest="summary_isp_txt",
        type=argparse.FileType("w"),
        default="_summary_isp.txt",
        help="Summary statistics on putative i-spanins",
    )
    parser.add_argument(
        "--putative_isp_gff",
        dest="putative_isp_gff",
        type=argparse.FileType("w"),
        default="_putative_isp.gff3",
        help="gff3 output for putative i-spanins",
    )

    parser.add_argument(
        "--max_isp",
        dest="max_isp",
        default=230,
        help="Maximum size of the ISP",
        type=int,
    )

    parser.add_argument(
        "--isp_mode",
        action="store_true",
        default=True
    )

    parser.add_argument(
        "--peri_min",
        type=int,
        default=18,
        help="amount of residues after TMD is found min"
    )

    parser.add_argument(
        "--peri_max",
        type=int,
        default=206,
        help="amount of residues after TMD is found max"
    )
    # parser.add_argument('-v', action='version', version='0.3.0') # Is this manually updated?
    args = parser.parse_args()
    the_args = vars(parser.parse_args())

    ### isp output, naive ORF finding:
    isps = OrfFinder(args.table, args.ftype, args.ends, args.isp_min_len, args.strand)
    isps.locate(
        args.fasta_file,
        args.out_isp_nuc,
        args.out_isp_prot,
        args.out_isp_bed,
        args.out_isp_gff3,
    )
    """
    >T7_EIS MLEFLRKLIPWVLVGMLFGLGWHLGSDSMDAKWKQEVHNEYVKRVEAAKSTQRAIGAVSAKYQEDLAALEGSTDRIISDLRSDNKRLRVRVKTTGISDGQCGFEPDGRAELDDRDAKRILAVTQKGDAWIRALQDTIRELQRK
    >lambda_EIS MSRVTAIISALVICIIVCLSWAVNHYRDNAITYKAQRDKNARELKLANAAITDMQMRQRDVAALDAKYTKELADAKAENDALRDDVAAGRRRLHIKAVCQSVREATTASGVDNAASPRLADTAERDYFTLRERLITMQKQLEGTQKYINEQCR
    """
    print(args.isp_mode)
    args.out_isp_prot.close()
    args.out_isp_prot = open(args.out_isp_prot.name, "r")

    pairs = tuple_fasta(fasta_file=args.out_isp_prot)

    # print(pairs)

    have_tmd = []  # empty candidates list to be passed through the user input criteria

    for (
        each_pair
    ) in (
        pairs
    ):  # grab transmembrane domains from spaninFuncts (queries for lysin snorkels # and a range of hydrophobic regions that could be TMDs)
        if len(each_pair[1]) <= args.max_isp:
            try:
                have_tmd += find_tmd(
                    pair=each_pair,
                    minimum=args.isp_min_dist,
                    maximum=args.isp_max_dist,
                    TMDmin=args.min_tmd_size,
                    TMDmax=args.max_tmd_size,
                    isp_mode=args.isp_mode,
                    peri_min=args.peri_min,
                    peri_max=args.peri_max,
                )
            except TypeError:
                continue

    if args.switch == "all":
        pass
    else:
        # for each_pair in have_lipo:
        range_of = args.switch
        range_of = re.search(("[\d]+:[\d]+"), range_of).group(0)
        start = int(range_of.split(":")[0])
        end = int(range_of.split(":")[1])
        have_tmd = parse_a_range(pair=have_tmd, start=start, end=end)

    total_isp = len(have_tmd)

    # ORF = [] # mightttttttttttt use eventually
    length = []  # grabbing length of the sequences
    candidate_dict = {k: v for k, v in have_tmd}
    with args.putative_isp_fa as f:
        for desc, s in candidate_dict.items():  # description / sequence
            f.write(">" + str(desc))
            f.write("\n" + lineWrapper(str(s).replace("*","")) + "\n")
            length.append(len(s))
            # ORF.append(desc)

    bot_size = min(length)
    top_size = max(length)
    avg = (sum(length)) / total_isp
    med = median(length)

    with args.summary_isp_txt as f:
        f.write("total potential o-spanins: " + str(total_isp) + "\n")
        f.write("average length (AA): " + str(avg) + "\n")
        f.write("median length (AA): " + str(med) + "\n")
        f.write("maximum orf in size (AA): " + str(top_size) + "\n")
        f.write("minimum orf in size (AA): " + str(bot_size))

    # Output the putative list in gff3 format
    args.putative_isp_fa = open(args.putative_isp_fa.name, "r")
    gff_data = prep_a_gff3(fa=args.putative_isp_fa, spanin_type="isp")
    write_gff3(data=gff_data, output=args.putative_isp_gff)

    """https://docs.python.org/3.4/library/subprocess.html"""
    """https://github.tamu.edu/CPT/Galaxy-Tools/blob/f0bf4a4b8e5124d4f3082d21b738dfaa8e1a3cf6/tools/phage/transmembrane.py"""
