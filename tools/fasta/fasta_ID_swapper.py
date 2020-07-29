
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

def swap_fa(fa):  # recycled from spaninFuncs.py
    """
        Splits a fasta into a tuple; (description, sequence)
        ---> INPUT : fasta file
        ---> OUTPUT : tuple of fasta separation
    """
    fasta = SeqIO.parse(fa, "fasta")
    descriptions = []
    sequences = []
    for r in fasta:  # iterates and stores each description and sequence
        description = r.description
        find_last_space = description.split(' ')[-1]
        remove_bracket = find_last_space.replace(']',' ')
        description = remove_bracket + description
        #print(find_last_space)
        #print(remove_bracket)
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
        "--swap_fa",
        type=argparse.FileType("w"),
        default="_id_swapped.fa",
        help="ID swapped fasta output",
    )  # appended file
    args = parser.parse_args()
    
    fasta = args.fasta
    new_headers = swap_fa(fasta)
    with args.swap_fa as f:
        for header in new_headers:
            f.write(">" + str(header[0]))  # fasta file 1 header_ fasta file 2 header
            f.write("\n" + line_wrapper(str(header[1])) + "\n")  # fasta file 1 seq + fasta file 2 seq