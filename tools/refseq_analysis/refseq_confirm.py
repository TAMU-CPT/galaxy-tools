import argparse
import re

#  CPT tools
from refseq_tools import CPTrefseq

def parse_inputs(text, file):
    """
    parse inputs (or combination of)
    """

    accs = []

    if text:
        if re.search(("__cn__"),str(text[0])):
            acc = text[0]
            split = acc.split("__cn__")
            accs.extend(split)
        else:
            accs.extend(text)
    
    if file:
        a = open(file.name).read().splitlines()
        accs.extend(a)
    
    if not accs:
        raise Exception("No accessions used in text box or file input, try again!")
    else:
        return accs

def parse_email(email):
    """
    Will likely become obsolete, but will join admin emails with user input email
    """

    ADMINS = ["curtisross@tamu.edu","cory.maughmer@tamu.edu","anthonyc@tamu.edu"]
    sep = ";"

    try:
        if "__at__" in email:
            split = email.split("__at__")
            email = f"{split[0]}@{split[1]}"
    except TypeError:
        raise Exception("Please Insert Email Address")

    ADMINS.insert(0,email)
    emails = sep.join(ADMINS)

    return emails

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Retrieves confirmation of existing RefSeq record with NCBI genome accession (gacc) or gid"
    )

    parser.add_argument(
        "--acc_file",
        type=argparse.FileType('r'),
        help="New line separated gacc and/or gid file",
    )

    parser.add_argument(
        "--acc_text",
        nargs="*",
        help="Accessions, separated by __cn__ for Galaxy, or space for offline use"
    )

    parser.add_argument(
        "--email",
        type=str,
        help="User email which will have CPT admins appended to list of emails for NCBI"
    )

    parser.add_argument(
        "--out_file",
        type=argparse.FileType('w'),
        default="_RefSeq_results.txt",
        help="Output tabular file"
    )
    args = parser.parse_args()

    #  Parse accessions from text of file
    accs = parse_inputs(args.acc_text, args.acc_file)

    #  Emails
    emails = parse_email(args.email)

    #  RefSeq confirmation pipeline execution
    report_list = []
    for acc in accs:
        og_acc = acc
        if acc: #  checks if not empty
            if type(acc) is str: 
                acc = acc.strip()
                acc = acc.split(".")[0]
            print(f"\nAccession: {acc}")
            depot = {
                "email" : emails,
                "acc" : acc
            }
            ref_acc, name = CPTrefseq(**depot).check_if_ref()
            temp_list = [og_acc,acc,ref_acc,name]
        else:
            temp_list = ["None","None","None","None"]
        report_list.append(temp_list)
    
    with args.out_file as file:
        file.writelines("Input_Accession\tQueried_Accession\tRefSeq_Match\tName\n")
        for each_return in report_list:
            file.writelines(f"{each_return[0]}\t{each_return[1]}\t{each_return[2]}\t{each_return[3]}\n")