import argparse
import os
import re
import time

from Bio import SeqIO

from excel_parser import ExcelParsing
import eutils


def is_binary_file(parser_var, arg):
    """
    THIS IS NOT WORKING AS DESIRED, MIGHT SCRAP AND FORCE NONBINARY EXCEL UPLOAD
    """
    if not os.path.exists(arg):
        parser_var.error(f"File {arg} does NOT exist!")
    else:
        try:
            return open(arg, 'rb')
        except:
            return open(arg, 'r')

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Prepares a excel file for a pairwise MIST dot plot")

    #  prep arguments
    parser.add_argument("excel_file", type=lambda x: is_binary_file(parser, x), help='Input excel file')
    parser.add_argument("--acc_col", type=str, help="column header label for accessions")
    parser.add_argument("--name_col", type=str, help="column header for MIST plot labels")
    parser.add_argument("--use_name_col", action="store_true", help="Uses column value for renaming the header")

    # eutils params
    parser.add_argument('--user_email', help="User email")
    parser.add_argument('--admin_email', help="Admin email")
    parser.add_argument('--history_file', help='Fetch results from previous query')
    parser.add_argument('--db', default="nuccore", help='Database to use')
    parser.add_argument('--api_key', help="NCBI API Key")
    #"65da3234a0dd70611ede507979d1f3885608"
    parser.add_argument('--retmode', default="fasta", help='Retmode')
    parser.add_argument('--rettype', default="fasta", help='Rettype')
    
    # output
    parser.add_argument('--output_fasta', type=argparse.FileType('w'), default="_MIST_multi.fa")

    args = parser.parse_args()

    #  parse data into dataframe using excel_parser
    cols = [str(args.acc_col).strip(), str(args.name_col).strip()]
    data = ExcelParsing(args.excel_file).chop_frame(cols=cols)
    
    #  prettify future headers
    names = list(data[args.name_col])
    spliced_names = []
    for name in names:
        if args.use_name_col: # just use what is in the column
            r = name.split(' ')
            spliced_names.append(r)
        else:
            s = name.split(' ')
            if re.search(('phage|virus|coli'), s[1]):
                r = s[2:]
            else:
                r = s[1:]
            spliced_names.append(r)
    ids = list(data[args.acc_col])
    combined_data = zip(spliced_names, ids)
    c = eutils.Client(
        history_file=args.history_file,
        user_email=args.user_email,
        admin_email=args.admin_email,
        api_key=args.api_key
    )

    #  retrieve data using accession column
    payload = {}
    for attr in ('retmode', 'rettype'):
        if getattr(args, attr, None) is not None:
            payload[attr] = getattr(args, attr)

    with args.output_fasta as f:
        for org in combined_data:
            payload['id'] = org[1]
            print(payload)
            obj = c.fetch(args.db, ftype=args.retmode, read_only_fasta=True, **payload)
            #print(obj)
            obj.description = obj.id
            obj.id = '_'.join(org[0])
            #print(obj.description)
            SeqIO.write(obj,"_temp.fa","fasta")
            for line in open("_temp.fa"):
                f.write(line)

