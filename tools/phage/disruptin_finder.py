'''
This program is intended to find genes with products that have a net charge greater than a given threshold
(default = +4) and length less than a given threshold (default = 80). This is a step in the process to find
genes similar to gene 28 in phage PhiKT - the first disruptin gene.

Input a multi fasta file with all of the predicted protein sequences from the genome as well as a threshold
size and charge. The program outputs a table and another fasta file. The table includes the name of the gene,
the net charge, length, total charged residues and the charge to size ratio for each product meeting the size
and charge criteria. The output fasta file includes records for all the sequences meeting the size and charge
criteria.

Example output:
Gene Name	Net Charge	Length	Number of Charged Residues	Charge to Size Ratio
F402_gp28	7	        56	    15	                        0.268

'''
import argparse
from Bio import SeqIO

def disruptin_finder(fasta_file, thresh_charge, thresh_size):
    # Creating output files - in tabular format
    print('Gene Name\tNet Charge\tLength\tNumber of Charged Residues\tCharge to Size Ratio')


    # Setting variable to calculate the net charge and total # of charged residues
    NetCharge = 0
    ChargeRes = 0

    total_record = []
    # iterates through each record within the fasta file
    for rec in SeqIO.parse(fasta_file,"fasta"):
        # Stores the name and sequence of each record
        Name, Sequence = rec.id, str(rec.seq)




        # For rec with length <= to the user-specified threshold, net charge and total charged residues are counted
        if len(Sequence) <= thresh_size:
            for aa in Sequence:
                # For R and K residues a positive charge is given
                if aa == 'R' or aa == 'K':
                    NetCharge += 1
                    ChargeRes += 1
                # For D and E residues a negative charge is given
                elif aa == 'D' or aa == 'E':
                    NetCharge -= 1
                    ChargeRes += 1
            # Charge to size ratio is calculated
            Length = len(Sequence)
            ChargeToSize = ChargeRes/Length

            # Prints the products that fit the criteria for a disruptin (size less than ~80 aa and net charge >+4) into
            # tabular format and outputs records to write a .fasta file
            if NetCharge >= thresh_charge:
                print('%s\t%d\t%d\t%d\t%0.3f' % (Name, NetCharge, Length, ChargeRes, ChargeToSize))
                total_record = total_record + [rec]


            # resetting the loop
            NetCharge = 0
            ChargeRes = 0

    yield total_record



if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    parser = argparse.ArgumentParser(description='Disruptin Finder')
    parser.add_argument('fasta_file', type=argparse.FileType("r"), help='Multi-FASTA Input')
    parser.add_argument('--thresh_charge', type=int, default=4)
    parser.add_argument('--thresh_size', type=int, default=80)
    args = parser.parse_args()

    for seq in disruptin_finder(**vars(args)):
        SeqIO.write(seq, sys.stdout, "fasta")
