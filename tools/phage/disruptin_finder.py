"""
This program is intended to find gene products that would be acceptable disruptin candidates.

The criteria can be toggled between selecting for proteins with:
    - net charge above a give threshold (default = +4) and length less than given threshold (default = 100 aa)
    OR
    - ratio of number of charged residues to length of the sequence above a given threshold (default = 0.25 residue/aa)
    and length less than given threshold (default = 100 aa)
    OR
    - net charge above a give threshold (default = +4), ratio of number of charged residues to length of the sequence
    above a given threshold (default = 0.25 residue/aa), and length less than given threshold (default = 100 aa)

Net charge of a sequence is calculated so that for every R or K residue the net charge increases by one, and for every
D or E residue the net charge decreases by one. The ratio of charged residues to length is calculated in a similar manner.
The residues R, K, D, and E each increase the number of charged residues by one, and total for the sequence is then
divided by the length to get the ratio.

Input a multi fasta file with all of the predicted protein sequences from the genome as well as a threshold
sequence length, net charge, and charge residue to length ratio. The program outputs another fasta file.
The output fasta file includes records for all the sequences meeting the size and charge criteria.

"""

from Bio import SeqIO
import argparse


def disruptin_finder(fasta_file, thresh_size, thresh_net_charge, thresh_charge_ratio, selection_criteria):
    # Iterable variables
    net_charge = 0
    charge_res = 0

    # Create record variable to store record information
    total_record = []

    # Parse the .fasta file and get the sequence
    for rec in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(rec.seq)

        if len(sequence) <= thresh_size:
            for aa in sequence:
                # For R and K residues a positive charge is given
                if aa in 'RK':
                    net_charge += 1
                    charge_res += 1
                # For D and E residues a negative charge is given
                elif aa in 'DE':
                    net_charge -= 1
                    charge_res += 1

            # Charge (total charged residues) to size ratio is calculated
            Length = len(sequence)
            charge_ratio = charge_res/Length

            # Based on the user-specified selection criteria a list of records is compiled
            if selection_criteria == 'net':
                if net_charge >= thresh_net_charge:
                    total_record = total_record + [rec]
            elif selection_criteria == 'ratio':
                if charge_ratio >= thresh_charge_ratio:
                    total_record = total_record + [rec]
            elif selection_criteria =='both':
                if charge_ratio >= thresh_charge_ratio and net_charge >= thresh_net_charge:
                    total_record = total_record + [rec]

            # Reset the iterable variables
            net_charge = 0
            charge_res = 0

    # The total list of records is returned by the function
    yield total_record


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    parser = argparse.ArgumentParser(description='Disruptin Finder')
    parser.add_argument('fasta_file', type=argparse.FileType("r"), help='Multi-FASTA Input')
    parser.add_argument('--thresh_net_charge', type=int, default=4)
    parser.add_argument('--thresh_size', type=int, default=100)
    parser.add_argument('--thresh_charge_ratio', type=float, default=0.25)
    parser.add_argument('--selection_criteria', action='store_true')
    args = parser.parse_args()

    for seq in disruptin_finder(**vars(args)):
        SeqIO.write(seq, sys.stdout, "fasta")
