import argparse

import pandas as pd
import re
from typing import List, Match, Dict, TextIO, Union
from datetime import date
# In this file, a "feature" refers to the collection of data between the > keys of the bprom output.
# That collection of data refers to one section of the DNA upstream of a gene

def read_bprom_file(bprom_file: TextIO) -> List[str]:
    """Reads in file, creating a list of strings with each list element containing a line from the file"""
    contents: List[str] = []

    with open(bprom_file) as file:
        for line in file:
            contents.append(line)

    return contents


def concatenate_then_split(contents: List[str]) -> List[str]:
    """Concatenates the file into one large string.
       Then splits it on '>' so that each feature's data is together in one element"""
    # Concatenates the entire file into one large string
    concat_contents: str = ''.join(contents)

    # Removing the empty string '' at element 0 used to make the join
    concat_contents: str = concat_contents[1:]

    # Splits the file into a list of strings on ">"
    features: List[str] = concat_contents.split('>')

    return features


def remove_promoterless_features(features: List[str]) -> List[str]:
    """For each concatenated feature string passed, removes the element
       if the # of predicted promoters is 0."""
    cleaned_features = features
    indices_to_delete = []
    for i, feature in enumerate(cleaned_features):
        if "Number of predicted promoters -      0" in cleaned_features[i]:
            indices_to_delete.append(i)
    # Must delete in reverse order, otherwise it changes the list indices after
    # the element deleted, and you delete subsequent elements at i+1, i+2, i+3, etc
    for i in sorted(indices_to_delete, reverse=True):
        del cleaned_features[i]

    return cleaned_features


def extract_accession(feature: str) -> str:
    """Extract accession"""
    accession: Match = re.search('[\w](.*)(?=_)', feature)
    accession: str = accession.group().replace('_', '').strip()

    return accession


def extract_test_seq_position(feature: str) -> List[str]:
    """Extract position in genome. Gets any number of values '(.*)' between the brackets
    using 'lookbehind/lookright' (?<=PATTERN) and 'lookahead/lookleft' regex assertions
    to extract (?<=Location=\\[)(.*)(?=]\\()"""
    location: Match = re.search('(?<=Location=\\[)(.*)(?=]\\()', feature)
    location: List[str] = location.group().split(':')

    return location


def extract_strand_direction(feature: str) -> str:
    """Extract strand direction for a feature, - or +"""
    # Matches for '(.)'
    direction: Match = re.search('(?<=\\().(?=\\))', feature)
    direction: str = direction.group()

    return direction


def extract_promoter_data(feature: str) -> Dict[str, str]:
    """Extracts all promoter data using regular expressions.
    Use for one element in the output of concatenate_then_split()"""
    # Extract promoter -10 and -35 sequences and scores
    # Gets everything between "-xx box at pos." and " Score"
    minus10: Match = re.search('(?<=-10 box at pos.)(.*)(?= Score)(.*)', feature)
    minus35: Match = re.search('(?<=-35 box at pos.)(.*)(?= Score)(.*)', feature)

    # Extracts the match and removes leading and trailing whitespace (which can be variable)
    # (the bprom output does not maintain the same # of whitespace characters
    # if there are less digits, at least for the scoring)
    minus10: List[str] = minus10.group().lstrip().split(' ')
    minus10_pos: int = int(minus10[0])
    minus10_seq: str = minus10[1]
    minus10_score: str = minus10[-1]

    minus35: List[str] = minus35.group().lstrip().split(' ')
    minus35_pos: int = int(minus35[0])
    minus35_seq: str = minus35[1]
    minus35_score: str = minus35[-1]

    # Can change these keys to change the column 9
    promoter_data: Dict[str, Union[str, int]] = {
        'minus10_pos': minus10_pos,
        'minus10_seq': minus10_seq,
        'minus10_score': minus10_score,
        'minus35_pos': minus35_pos,
        'minus35_seq': minus35_seq,
        'minus35_score': minus35_score
    }
    return promoter_data


def convert_extracted_promoter_data_to_ID_column_format(
        promoter_data: Dict[str, Union[str, int]],
        calculated_promoter_positions: List[int]) -> str:
    """Converts input data to the GFF3 ID column (column 9) format, a semicolon separated
       list of values providing additional information about each feature"""
    # Replaces the BPROM output positions with the calculated ones
    minus_10_calculated: int = calculated_promoter_positions[2]
    minus_35_calculated: int = calculated_promoter_positions[3]
    promoter_data['minus10_pos'] = minus_10_calculated
    promoter_data['minus35_pos'] = minus_35_calculated

    # Creates the column 9 string (attributes)
    promoter_data: List[Union[str, int]] = [f'{key}={value}' for key, value in promoter_data.items()]
    promoter_data: str = 'Description=Predicted promoter data;' + 'Note=' + ','.join(promoter_data) + ';'
    return promoter_data


def extract_LDF_score(feature: str) -> str:
    """Extract LDF score"""
    LDF: Match = re.search('(?<=LDF-)(.*)', feature)
    LDF: str = LDF.group().strip()

    return LDF


def calculate_promoter_position(feature: str):
    """Calculate promoter positions (in the context of the genome) based on BPROM predictions."""
    # Get 'Promoter Pos:     X' data. This refers to the predicted transcriptional start site!
    promoter_pos: Match = re.search('(?<=Promoter Pos:)(.*)(?=LDF)', feature)
    promoter_pos: int = int(promoter_pos.group().strip())

    # Get start and end positions from 'Location=[XXX:YYYY]'
    test_seq_position: List[str] = extract_test_seq_position(feature)
    test_cds_location_start_pos: int = int(test_seq_position[0])
    test_cds_location_end_pos: int = int(test_seq_position[1])

    promoter_data: Dict[str, Union[str, int]] = extract_promoter_data(feature)

    ''' IMPORTANT!! Whether or not you add or subtract to calculate the promoter start
    # position depends on whether we're on the + or - strand!
    # The workflow Jolene uses is smart enough to correctly pull upstream
    # for both + and - strands (i.e., pulls left for +, pulls right for -)
    # THEREFORE, for a gene with a start at 930 on the + strand, it pulls 830:930
    # And for a gene with a start at 930 on the - strand, it pulls 930:1030 '''

    direction: str = extract_strand_direction(feature)

    if direction == '+':
        # BPROM starts counting from the LEFT boundary for + strand test sequences (as expected)
        # Get -10 promoter position
        minus10_pos: int = promoter_data['minus10_pos']
        minus10_pos_in_context_of_genome: int = test_cds_location_start_pos + minus10_pos
        # Get -35 promoter position
        minus35_pos: int = promoter_data['minus35_pos']
        minus35_pos_in_context_of_genome: int = test_cds_location_start_pos + minus35_pos

        start: int = test_cds_location_start_pos + minus35_pos
        end: int = test_cds_location_start_pos + promoter_pos

        calculated_promoter_positions: List[int] = [
            start, end, minus10_pos_in_context_of_genome, minus35_pos_in_context_of_genome]
        return calculated_promoter_positions

    elif direction == '-':
        # BPROM starts counting from the RIGHT boundary for - strand test sequences
        # Get -10 promoter position
        minus10_pos: int = promoter_data['minus10_pos']
        minus10_pos_in_context_of_genome: int = test_cds_location_end_pos - minus10_pos
        # Get -35 promoter position
        minus35_pos: int = promoter_data['minus35_pos']
        minus35_pos_in_context_of_genome: int = test_cds_location_end_pos - minus35_pos

        # The start and end are reversed
        end: int = test_cds_location_end_pos - minus35_pos
        start: int = test_cds_location_end_pos - promoter_pos

        calculated_promoter_positions: List[int] = [
            start, end, minus10_pos_in_context_of_genome, minus35_pos_in_context_of_genome]
        return calculated_promoter_positions

    else:
        assert "Error: Strand data neither \'+\' nor \'-\'"


def extract_tf_binding_elements():
    """Extract predicted transcription factor binding elements"""
    return


def extract_data_for_all_features(features: List[str]) -> List[List[Union[str, int]]]:
    """Loops through cleaned bprom output extracting all data of interest and builds the
       structure for loading into a dataframe"""
    extracted_data: List[List[Union[str, int]]] = []
    for feature in features:
        # loop through features, a List[str] containing each feature [str] in the
        # original bprom format as a single string, but cleaned of irrelevant data
        calculated_promoter_positions: List[int] = calculate_promoter_position(feature)
        promoter_data: Dict[str, str] = extract_promoter_data(feature)
        promoter_data_converted: str = convert_extracted_promoter_data_to_ID_column_format(
            promoter_data, calculated_promoter_positions)

        extracted_data.append(
            [extract_accession(feature),  # Seqid, col 1
             'bprom',  # Source, col 2
             'promoter',  # Type, col 3
             calculated_promoter_positions[0],  # Start, col 4
             calculated_promoter_positions[1],  # End, col 5
             extract_LDF_score(feature),  # Score, col 6
             extract_strand_direction(feature),  # Strand direction, col 7
             '.',  # Phase, col 8
             promoter_data_converted,  # Attributes, col 9
             ])

    return extracted_data


def convert_to_dataframe(extracted_data: List[List[str]]) -> pd.DataFrame:
    """Convert extracted and processed data to Pandas dataframe with gff3 column names"""

    df = pd.DataFrame(extracted_data,
                      columns=['seqid', 'source', 'type', 'start', 'end',
                               'score', 'strand', 'phase', 'attributes']
                      )
    return df


def write_to_gff3(dataframe) -> None:
    """Create a gff3 text file from the DataFrame by converting to a tab separated values (tsv) file"""
    tsv: pd.DataFrame = dataframe.to_csv(sep='\t', index=False, header=None)

    # Gets the first element of the first column to use for
    accession: str = dataframe.iloc[0][0]

    year, month, day = date.today().year, date.today().month, date.today().day

    #with open(f'{year}_{month}_{day}_bprom_as_gff3_{accession}.txt', 'w') as wf:
        # Header so Galaxy can recognize as GFF3
    print('##gff-version 3\n')
    for line in tsv:
      print(line)

    return

def convert_bprom_output_to_gff3(bprom_file: TextIO) -> None:
    """Master function. Given a BPROM .txt file as output, extracts data and writes as a GFF3 file"""
    bprom_file: List[str] = read_bprom_file(bprom_file)
    concatenated_bprom_file: List[str] = concatenate_then_split(bprom_file)
    working_file: List[str] = remove_promoterless_features(concatenated_bprom_file)
    extracted_data: List[List[Union[str, int]]] = extract_data_for_all_features(working_file)
    gff3_dataframe: pd.DataFrame = convert_to_dataframe(extracted_data)
    # Create the gff3 text file
    write_to_gff3(gff3_dataframe)

    return


if __name__ == '__main__':
    ## Shows the DataFrame output in the terminal for testing/debugging
    # bprom_file = read_bprom_file('BPROM_output.txt')
    # concatenated_bprom_file: List[str] = concatenate_then_split(bprom_file)
    # working_file = remove_promoterless_features(concatenated_bprom_file)
    # print(convert_to_dataframe(extract_data_for_all_features(working_file)).to_string())

    parser = argparse.ArgumentParser(
        description='converts BPROM output to the gff3 file format')

    parser.add_argument('-f', help='bprom file as .txt')
    args = parser.parse_args()
    # Actual function for converting the BPROM output to gff3
    convert_bprom_output_to_gff3(args.f)

    # Upload to cpt github in the directory Galaxy-Tools/tools/external/
