import argparse
import os
import re

from excel_parser import ExcelParsing


def is_binary_file(parser, arg):
    if not os.path.exists(arg):
        parser.error(f"File {arg} does NOT exist!")
    else:
        try:
            return open(arg, 'rb')
        except:
            return open(arg, 'r')

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Prepares a excel file for a pairwise MIST dot plot")

    parser.add_argument("excel_file", type=lambda x: is_binary_file(parser, x), help='Input excel file')
    parser.add_argument("--acc_col", type=str, help="column header label for accessions")
    parser.add_argument("--name_col", type=str, help="column header for MIST plot labels")
    parser.add_argument("--output_multifa", type=argparse.FileType('w'), default="_multi.fa")

    #parser.add_argument("--cols", type=)

    args = parser.parse_args()

    #  parse data into dataframe using excel_parser
    cols = [args.acc_col, args.name_col]
    data = ExcelParsing(args.excel_file).chop_frame(cols=cols)

    #  prettify future headers
    names = list(data[args.name_col])
    spliced_names = []
    for name in names:
        s = name.split(' ')
        if re.search(('phage|virus|coli'), s[1]):
            r = s[2:]
        else:
            r = s[1:]
        spliced_names.append(r)
    
    print(spliced_names)
    #  retrieve data using accession column
    """
    accs = list(data[args.acc_col])
    for acc in accs:
        load = {
            "email":"curtisross@tamu.edu",
            "acc":acc,
            "db":"nuccore",
            "ret_type":"fasta"
        }
        c = SeqIO.read(CPTEfetch(**load).retrieve_data(), "fasta")
        #c = str(c)
        print(c)
        exit()
    """
