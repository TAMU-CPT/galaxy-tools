
from Bio import SeqIO
from Bio import Seq
import argparse

#  Takes a multifasta file and appends the strain (inside closing brackets) to front
#  Use is for MSA

def line_wrapper(text, charactersize=60):
    """
    Used to wrap sequences for fasta output
    """

    if len(text) <= charactersize:
        return text
    else:
        return (
            text[:charactersize]
            + "\n"
            + line_wrapper(text[charactersize:], charactersize)
        )

def swap_fa(fa, wrapper="bracket", space="last"):  # recycled from spaninFuncs.py
    """
    Splits a fasta into a tuple; (description, sequence)
    ---> INPUT : fasta file
    ---> OUTPUT : tuple of fasta separation
    space is the specific spot a user wants to be the new leading word in the header
    wrapper is either "bracket" or "parenth" for different vars
    """
    if space != "last": #  just bringing everything in as str for simplification from Galaxy...
        space = int(space)
    fasta = SeqIO.parse(fa, "fasta")
    descriptions = []
    sequences = []
    if wrapper == "bracket":
        wrap = ['[',']']
    elif wrapper == "parenth": #  can expand more if desired
        wrap = ['(',')']
    for r in fasta:  # iterates and stores each description and sequence
        description = r.description
        if space == "last":
            find_space = description.split(' ')[-1]
            remove_wrap = find_space.replace(str(wrap[1]),' ')
        else:
            if description.split(' ')[-1] == description.split(' ')[space-1]: # means we're at the end of the header
                find_space = description.split(' ')[-1]
                remove_wrap = find_space.replace(str(wrap[1]),' ')
            else:
                find_space = description.split(' ')[space-1]
                if find_space[-1] == wrap[0]:
                    remove_wrap = find_space.replace(str(wrap[0]),' ')
                else:
                    remove_wrap = find_space + ' '
        description = remove_wrap + description
        sequence = str(r.seq)
        descriptions.append(description)
        sequences.append(sequence)
    return zip(descriptions, sequences)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="append two fasta files into one")
    parser.add_argument(
        "fasta", type=argparse.FileType("r"), help="input multifasta file"
    )  # the first "input" argument
    parser.add_argument(
        "--wrapper",
        default="bracket",
        choices=("bracket","parenth"),
        help="Ending (or starting) bracket or parentheses to remove."
    )
    parser.add_argument(
        "--space",
        default="last",
        type=str,
        help="Space number to append to front of header"
    )
    parser.add_argument(
        "--swap_fa",
        type=argparse.FileType("w"),
        default="_id_swapped.fa",
        help="ID swapped fasta output",
    )  # appended file
    args = parser.parse_args()
    
    fasta = args.fasta
    new_headers = swap_fa(fasta, wrapper=args.wrapper, space=args.space)
    with args.swap_fa as f:
        for header in new_headers:
            f.write(">" + str(header[0]))  # fasta file 1 header_ fasta file 2 header
            f.write("\n" + line_wrapper(str(header[1])) + "\n")  # fasta file 1 seq + fasta file 2 seq