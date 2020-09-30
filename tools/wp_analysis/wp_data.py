import argparse
import re

#  CPT tools
from wp_tools import CPTLink

def parse_inputs(text,file,galaxy_mode=False):
    """
    Parses the inputs of a text box and new line separated pacc file
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
        raise Exception("No accessions used, check file and input.")
    else:
        return accs

def parse_email(email):
    """
    Parses user input email and appends CPT Admins to NCBI email
    """
    ADMINS = ["curtisross@tamu.edu","cory.maughmer@tamu.edu","anthonyc@tamu.edu"]
    sep = ';'

    try:
        if "__at__" in email:
            split = email.split("__at__")
            email = f"{split[0]}@{split[1]}"
    except TypeError:
        raise Exception("Please Insert Email Address")

    ADMINS.insert(0,email)
    emails = sep.join(ADMINS)

    return emails

def write_table(list_of_data, file):
    """
    writes output table, uses list of data from CPTlink output
    """
    with file as f:
        f.write(f"WP_accession\tGenome_accession\tTaxID\tOrganism\tWP_count\n")
        for acc in list_of_data:
            for gacc_data in acc[1]:
                f.write(f"{acc[0]}\t{gacc_data[0]}\t{gacc_data[1]}\t{gacc_data[2]}\t{acc[2]}\n")

def write_gaccs(list_of_data, file):
    """
    writes output gacc file, uses list of data from CPTlink output
    """

    for acc in list_of_data:
        print(acc)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Retrieve information from a WP accession"
    )

    parser.add_argument(
        "--wp_text",
        nargs="*",
        help="WP accessions, separated by __cn__ for Galaxy, or space for offline"
    )

    parser.add_argument(
        "--wp_file",
        type=argparse.FileType("r"),
        help="New line separated WP accessions file"
    )

    parser.add_argument(
        "--email",
        type=str,
        help="Entrez requires an email to connect to NCBI database. CPT Admin emails will be appended to list."
    )

    parser.add_argument(
        "--wp_amt",
        dest='wp_amt',
        choices=('first','all'),
        default='first',
    )

    parser.add_argument(
        "--out_table",
        type=argparse.FileType("w"),
        default="_return_table.txt",
        help="Output table consisting of accession data"
    )

    args = parser.parse_args()

    #  Get accessions from input file and/or text
    accs = parse_inputs(args.wp_text,args.wp_file)
    
    #  Emails
    emails = parse_email(args.email)

    #  Run functions
    package = {
        "email" : emails,
        "db" : "nuccore",
        "dbfrom" : "protein",
    }
    wps = []
    for acc in accs:
        package["acc"] = acc
        if args.wp_amt == 'all': #  sorta a hacky way to check and see if we're grabbing first or all
            wp_all = True
        else:
            wp_all = False

        pacc, gacc, wp_amt = CPTLink(**package).map_accessions(wp_all)
        
        current_wp = [pacc,gacc,wp_amt]
        wps.append(current_wp)
    
    write_table(wps,args.out_table)

