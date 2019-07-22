#!/usr/bin/env python

import csv
import argparse
#import sys

# This function reads through the tagseq file and outputs a list of sequence names and the lengths of each sequence.
def garnier_sequences(tagseq_file = None):
    # open the file and create blank lists
    f = tagseq_file #open(tagseq_file, 'r')
    f.seek(0)
    sequence = []
    lengths = []

    # for each line the in file, search for the words 'Sequence' and 'to' to find the sequence name and length,
    # respectively. Then add sequence names and lengths to the proper lists
    for line in f:
        words = line.split()
        if line.startswith('# Sequence:'):
        #if 'Sequence:' in line:
            #if words[1] == 'Sequence:':
            sequence += [words[words.index('Sequence:') + 1]]
            #if words[5] == 'to:':
            #    lengths += [int(words[6])]
            if words.index('to:'):
                lengths += [int(words[words.index('to:') + 1])]
    # return the sequence names and lengths
    return sequence, lengths


# This function extracts the helix, sheet, turn, and coil predictions from the file. The predictions for each type of
# secondary structure are joined together in one string.
def garnier_secondary_struct(tagseq_file = None):
    # opens the file and sets variables for the structural predictions
    f = tagseq_file #open(tagseq_file, 'r')
    helix = ''
    turns = ''
    coil = ''
    sheet = ''

    # if the first work in the line indicates a structural prediction, it adds the rest of the line to the right
    # prediction string.
    for line in f:
        words = line.split()
        if len(words) > 0:
            if words[0] in 'helix':
                helix += str(line[6:]).rstrip('\n')
            elif words[0] in 'sheet':
                sheet += str(line[6:]).rstrip('\n')
            elif words[0] in 'turns':
                turns += str(line[6:]).rstrip('\n')
            elif words[0] in 'coil':
                coil += str(line[6:]).rstrip('\n')
    # f.close()
    # returns the four structural prediction strings
    return helix, turns, coil, sheet


# This functions cuts the strings based on the lengths of the original sequences. Lengths are given in a list.
def vector_cutter(vector, lengths_to_cut):
    # sets up iteration variables
    start = 0
    end = lengths_to_cut[0]
    maximum = len(lengths_to_cut)
    # creates output list
    output = []
    # loops through the number of sequences based on the number of lengths
    for i in range(maximum):
        # outputs list of sequence strings
        output += [str(vector[start:end])]
        start = end
        if i + 1 != maximum:
            end += lengths_to_cut[i + 1]
    # returns list of strings. Each sequence has a string included in the list.
    return output


# this function takes the helix, turn, sheet, and coil predictions for each sequence and creates a single structural
# prediction string.
def single_prediction(helix, sheet, turns, coil):
    # sets output list
    secondary_structure = []
    # checks to make sure each of the strings is the same length
    if len(helix) == len(sheet) == len(coil) == len(turns):
        # loops through the length of each sequence, and when the value is not a blank it is added to the output
        # prediction list.
        for j in range(len(helix)):
            if helix[j] != ' ':
                secondary_structure += [str(helix[j])]
            elif sheet[j] != ' ':
                secondary_structure += [str(sheet[j])]
            elif coil[j] != ' ':
                secondary_structure += [str(coil[j])]
            else:
                secondary_structure += [str(turns[j])]
    # returns the output prediction list for the sequence
    return secondary_structure



if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    parser = argparse.ArgumentParser(description='Read Garnier Secondary Structure Prediction')
    parser.add_argument('tagseq_file', type=argparse.FileType("r"), help='Tagseq file input')
    args = parser.parse_args()

    # opens the tagseq file and prepares for writing csv
    #f = open(sys.stdout, 'w', newline='')
    # writer = csv.writer(f)

    # reads tagseq file for helix, turn, coil, and sheet sequences as well as for names and lengths of the sequences
    # summarized in the tagseq file#!/usr/bin/env python\r
    Hel, Tur, Coi, She = garnier_secondary_struct(**vars(args))
    names, gives = garnier_sequences(**vars(args))

    # cut each of the structural prediction strings so that they are individual sequences
    Helix = vector_cutter(Hel, gives)
    Sheet = vector_cutter(She, gives)
    Turns = vector_cutter(Tur, gives)
    Coil = vector_cutter(Coi, gives)

    # for each sequence compile the four types of structural predictions into a single prediction, and output the final
    # prediction in csv format and to the screen
    for i in range(len(Helix)):
        Final = single_prediction(Helix[i], Sheet[i], Turns[i], Coil[i])
        #csv.writerow(['Sequence: '] + [names[i]])
        #csv.writerow(Final)
        print('Sequence Name: ' + '\t' + names[i])
        print('\t'.join(Final))
