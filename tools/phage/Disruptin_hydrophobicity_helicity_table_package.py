"""
This program is intended to create the output table for the disruptin finder workflow
"""
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import ProtParamData
import csv
import argparse


def disruptin_table(garnier_file, fasta_file):
    # Iterable variables
    position = 1
    net_charge = 0
    charge_res = 0
    record_number = 0

    # loop structures
    names = []
    sec_struct = []

    # reading the lines from the garnier csv file
    with open(garnier_file, "r") as csvfile:
        garnierreader = csv.reader(csvfile)
        for row in garnierreader:
            if row[0] == "Sequence: ":
                names += [row[1]]
            elif row[0] in "HETC":
                sec_struct += ["".join(row)]
    record = []
    p = []
    r = []
    c = []
    h = []
    s = []

    # Parse the .fasta file and get the sequence
    for rec in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(rec.seq)

        # Set up the information vectors: for position #, residue, hydrophobic/charge/polar/nonpolar, and secondary
        # structure
        record += [rec.id]
        position_vec = []
        residue_vec = []
        charge_sym_vec = []
        sec_struct_vec = []

        for aa in sequence:
            position_vec += [str(position)]
            residue_vec += [str(aa)]
            sec_struct_vec += [str(sec_struct[record_number][position - 1])]

            # For R and K residues a positive charge is given
            if aa in "RK":
                symbol = "+"
            # For D and E residues a negative charge is given
            elif aa in "DE":
                symbol = "-"
            elif aa in "AVMILPWFG":
                symbol = "N"
            elif aa in "HSYTCQN":
                symbol = "P"
            charge_sym_vec += symbol
            position += 1

            # Calculating hyrophobicity based on Kyte and Doolittle scale. Using binning value of 9. Since the binning
            # is 9, the first 4 residues and last 4 residues as set blank so as to center the values to their
            # approximate position on the sequence.
            prot_ana_seq = ProteinAnalysis(sequence)
            hydro = [0] * 4 + prot_ana_seq.protein_scale(ProtParamData.kd, 9) + [0] * 4

        record_number += 1
        position = 1

        p += [position_vec]
        r += [residue_vec]
        c += [charge_sym_vec]
        h += [hydro]
        s += [sec_struct_vec]

    # returns values for name of the sequence
    return record, p, r, c, h, s


if __name__ == "__main__":
    # Grab all of the filters from our plugin loader
    parser = argparse.ArgumentParser(description="Disruptin Table Output")
    parser.add_argument(
        "garnier_file", type=argparse.FileType("r"), help="csv file from garnier reader"
    )
    parser.add_argument(
        "fasta_file",
        type=argparse.FileType("r"),
        help="fasta file of disruptin candidates",
    )
    args = parser.parse_args()

    # Set up output location
    f = open(sys.stdout, "w", newline="")
    writer1 = csv.writer(f)

    iden, position, residue, charge, hydro, struct = disruptin_table(**vars(args))

    for i in range(len(iden)):
        writer1.writerow(["Protein ID"] + [iden[i]])
        writer1.writerow(["Position"] + [format(x, "s") for x in position[i]])
        writer1.writerow(["Residue"] + [format(x, "s") for x in residue[i]])
        writer1.writerow(["Charge"] + [format(x, "s") for x in charge[i]])
        writer1.writerow(["Hydrophobicity"] + [format(x, ".3f") for x in hydro[i]])
        writer1.writerow(["Secondary Structure"] + [format(x, "s") for x in struct[i]])
        writer1.writerow([""])

        print(str(iden[i]))
        print("Position \t " + "\t".join(position[i]))
        print("Residue \t" + "\t".join(residue[i]))
        print("Charge \t" + "\t".join(charge[i]))
        print("Hydrophobicity \t" + "\t".join(format(x, ".3f") for x in hydro[i]))
        print("Secondary Structure \t" + "\t".join(struct[i]))
