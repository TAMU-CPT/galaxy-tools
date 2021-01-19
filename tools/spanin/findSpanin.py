##### findSpanin.pl --> findSpanin.py
######### Much of this code is very "blocked", in the sense that one thing happens...then a function happens on the return...then another function...etc...etc...

import argparse
import os
import re # new
import itertools # new
from collections import Counter, OrderedDict
from spaninFuncs import getDescriptions, grabLocs, spaninProximity, splitStrands, tuple_fasta, lineWrapper

### Requirement Inputs
#### INPUT : putative_isp.fa & putative_osp.fa (in that order)
#### PARAMETERS :

###############################################################################
def write_output(candidates):
    """ output file function...maybe not needed """
    pass

def reconfigure_dict(spanins):
    """
    re organizes dictionary to be more friendly for checks
    """

    new_spanin_dict = {}

    for each_spanin_type, data_dict in spanins.items():
        #print(f"{each_spanin_type} == {data_dict}")
        new_spanin_dict[each_spanin_type] = {}
        new_spanin_dict[each_spanin_type]['positive'] = {}
        new_spanin_dict[each_spanin_type]['negative'] = {}
        new_spanin_dict[each_spanin_type]['positive']['coords'] = []
        new_spanin_dict[each_spanin_type]['negative']['coords'] = []
        for outter_orf, inner_data in data_dict.items():
            list_of_hits = []
            for data_content in inner_data:
                #print(data_content)
                data_content.insert(0, outter_orf)
                #print(f"new data_content -> {data_content}")
                #print(data_content)
                #list_of_hits += [data_content]
                #new_spanin_dict[each_spanin_type] += [data_content]
                if data_content[6] == "+":
                    #print(f"{each_spanin_type} @ POSITIVE")
                    new_spanin_dict[each_spanin_type]['positive']['coords'] += [data_content]
                elif data_content[6] == "-":
                    #print(f"{each_spanin_type} @ NEGATIVE")
                    new_spanin_dict[each_spanin_type]['negative']['coords'] += [data_content]
        #print(new_spanin_dict[each_spanin_type])
            #print(reorganized)
            #print(f"{outter_orf} => {inner_data}")
            #print(new_spanin_dict)
        
        #print('\n')
    #for k, v in new_spanin_dict.items():
        #print(k)
        #print(v)
    return new_spanin_dict



def check_for_uniques(spanins):
    """
    Checks for unique spanins based on spanin_type.
    If the positive strand end site is _the same_ for a i-spanin, we would group that as "1".
    i.e. if ORF1, ORF2, and ORF3 all ended with location 4231, they would not be unique.
    """
    pair_dict = {}
    pair_dict = {
        'pairs' : {
            'location_amount' : [],
            'pair_number' : {},
        }
    }
    for each_spanin_type, spanin_data in spanins.items():
        #print(f"{each_spanin_type} ===> {spanin_data}")
        # early declarations for cases of no results
        pos_check = [] # end checks
        pos_uniques = []
        neg_check = [] # start checks
        neg_uniques = []
        unique_ends = []
        pos_amt_unique = 0
        neg_amt_unique = 0
        amt_positive = 0
        amt_negative = 0
        spanin_data['uniques'] = 0
        spanin_data['amount'] = 0
        #spanin_data['positive']['amt_positive'] = 0
        #spanin_data['positive']['pos_amt_unique'] = 0
        #spanin_data['positive']['isp_match'] = []
        #spanin_data['negative']['amt_negative'] = 0
        #spanin_data['negative']['neg_amt_unique'] = 0
        #spanin_data['negative']['isp_match'] = []
        #print(spanin_data)
        if spanin_data['positive']['coords']:
            # do something...
            #print('in other function')
            #print(spanin_data['positive']['coords'])
            for each_hit in spanin_data['positive']['coords']:
                pos_check.append(each_hit[2])
                pair_dict['pairs']['location_amount'].append(each_hit[2])
                pos_uniques = list(set([end_site for end_site in pos_check if pos_check.count(end_site) >= 1]))
            #print(pos_check)
            #print(pos_uniques)
            amt_positive = len(spanin_data['positive']['coords'])
            pos_amt_unique = len(pos_uniques)
            if amt_positive:
                spanin_data['positive']['amt_positive'] = amt_positive
                spanin_data['positive']['pos_amt_unique'] = pos_amt_unique
            #pair_dict['pairs']['locations'].extend(pos_uniques)
        else:
            spanin_data['positive']['amt_positive'] = 0
            spanin_data['positive']['pos_amt_unique'] = 0
        if spanin_data['negative']['coords']:

            # do something else...
            #print('in other function')
            #print(spanin_data['negative']['coords'])
            for each_hit in spanin_data['negative']['coords']:
                neg_check.append(each_hit[1])
                pair_dict['pairs']['location_amount'].append(each_hit[1])
                neg_uniques = list(set([start_site for start_site in neg_check if neg_check.count(start_site) >= 1]))
            #print(neg_uniques)
            amt_negative = len(spanin_data['negative']['coords'])
            neg_amt_unique = len(neg_uniques)
            if amt_negative:
                spanin_data['negative']['amt_negative'] = amt_negative
                spanin_data['negative']['neg_amt_unique'] = neg_amt_unique
            #pair_dict['pairs']['locations'].extend(neg_uniques)
        else:
            spanin_data['negative']['amt_negative'] = 0
            spanin_data['negative']['neg_amt_unique'] = 0
        spanin_data['uniques'] += (spanin_data['positive']['pos_amt_unique'] + spanin_data['negative']['neg_amt_unique'])
        spanin_data['amount'] += (spanin_data['positive']['amt_positive'] + spanin_data['negative']['amt_negative'])
        #print(spanin_data['uniques'])
    list(set(pair_dict['pairs']['location_amount']))
    pair_dict['pairs']['location_amount'] = dict(Counter(pair_dict['pairs']['location_amount']))
    for data in pair_dict.values():
        #print(data['locations'])
        #print(type(data['locations']))
        v = 0
        for loc, count in data['location_amount'].items():
            #data['pair_number'] = {loc
            v += 1
            data['pair_number'][loc] = v
    #print(dict(Counter(pair_dict['pairs']['locations'])))
    #print(pair_dict)
    spanins['total_amount'] = spanins['EMBEDDED']['amount'] + spanins['SEPARATED']['amount'] + spanins['OVERLAPPED']['amount']
    #spanins['total_unique'] = spanins['EMBEDDED']['uniques'] + spanins['SEPARATED']['uniques'] + spanins['OVERLAPPED']['uniques']
    spanins['total_unique'] = len(pair_dict['pairs']['pair_number'])
    return spanins, pair_dict

if __name__ == "__main__":

    # Common parameters for both ISP / OSP portion of script

    parser = argparse.ArgumentParser(
        description="Trim the putative protein candidates and find potential i-spanin / o-spanin pairs"
    )

    parser.add_argument(
        "putative_isp_fasta_file",
        type=argparse.FileType("r"),
        help='Putative i-spanin FASTA file, output of "generate-putative-isp"',
    )  # the "input" argument

    parser.add_argument(
        "putative_osp_fasta_file",
        type=argparse.FileType("r"),
        help='Putative o-spanin FASTA file, output of "generate-putative-osp"',
    )

    parser.add_argument(
        "--max_isp_osp_distance",
        dest="max_isp_osp_distance",
        default=10,
        type=int,
        help="max distance from end of i-spanin to start of o-spanin, measured in AAs",
    )

    parser.add_argument(
        "--embedded_txt",
        dest="embedded_txt",
        type=argparse.FileType("w"),
        default="_findSpanin_embedded_results.txt",
        help="Results of potential embedded spanins",
    )
    parser.add_argument(
        "--overlap_txt",
        dest="overlap_txt",
        type=argparse.FileType("w"),
        default="_findSpanin_overlap_results.txt",
        help="Results of potential overlapping spanins",
    )
    parser.add_argument(
        "--separate_txt",
        dest="separate_txt",
        type=argparse.FileType("w"),
        default="_findSpanin_separated_results.txt",
        help="Results of potential separated spanins",
    )

    parser.add_argument(
        "--summary_txt",
        dest="summary_txt",
        type=argparse.FileType("w"),
        default="_findSpanin_summary.txt",
        help="Results of potential spanin pairs",
    )
    parser.add_argument(
        "-v", action="version", version="0.3.0"
    )  # Is this manually updated?
    args = parser.parse_args()


    #### RE-WRITE
    SPANIN_TYPES = {}
    SPANIN_TYPES['EMBEDDED'] = {}
    SPANIN_TYPES['OVERLAPPED'] = {}
    SPANIN_TYPES['SEPARATED'] = {}
    #SPANIN_TYPES = {
    #    'EMBEDDED' : {},
    #    'OVERLAPPED' : {},
    #    'SEPARATED' : {},
    #}

    isp = getDescriptions(args.putative_isp_fasta_file)
    args.putative_isp_fasta_file = open(args.putative_isp_fasta_file.name, "r")
    isp_full = tuple_fasta(args.putative_isp_fasta_file)

    osp = getDescriptions(args.putative_osp_fasta_file)
    args.putative_osp_fasta_file = open(args.putative_osp_fasta_file.name, "r")
    osp_full = tuple_fasta(args.putative_osp_fasta_file)

    #### location data
    location_data = {
        'isp' : [],
        'osp' : []
    }
    spanins = [isp, osp]
    for idx, each_spanin_type in enumerate(spanins):
        for description in each_spanin_type:
            locations = grabLocs(description)
            if idx == 0: # i-spanin
                location_data['isp'].append(locations)
            elif idx == 1: # o-spanin
                location_data['osp'].append(locations)

    #### Check for types of spanins
    embedded, overlap, separate = spaninProximity(
                                        isp=location_data['isp'],
                                        osp=location_data['osp'], 
                                        max_dist=args.max_isp_osp_distance * 3
                                    )
    
    SPANIN_TYPES['EMBEDDED'] = embedded
    SPANIN_TYPES['OVERLAPPED'] = overlap
    SPANIN_TYPES['SEPARATED'] = separate

    #for spanin_type, spanin in SPANIN_TYPES.items():
    #    s = 0
    #    for sequence in spanin.values():
    #        s += len(sequence)
    #    SPANIN_TYPES[spanin_type]['amount'] = s
    #    SPANIN_TYPES[spanin_type]['unique'] = len(spanin.keys())
    
    #check_for_unique_spanins(SPANIN_TYPES)
    spanins = reconfigure_dict(SPANIN_TYPES)
    spanins, pair_dict = check_for_uniques(spanins)
    #print(pair_dict)
    with args.summary_txt as f:
        for each_spanin_type, spanin_data in spanins.items():
            try:
                if each_spanin_type not in ["total_amount","total_unique"]:
                    #print(each_spanin_type)
                    #print(each_spanin_type)
                    f.write("=~~~~~= "+str(each_spanin_type) +" Spanin Candidate Statistics =~~~~~=\n")
                    f.writelines("Total Candidate Pairs = "+str(spanin_data['amount'])+"\n")
                    f.writelines("Total Unique Pairs = "+str(spanin_data['uniques'])+"\n")
                    if each_spanin_type == "EMBEDDED":
                        for k, v in SPANIN_TYPES['EMBEDDED'].items():
                            #print(k)
                            f.writelines(""+str(k)+" ==> Amount of corresponding candidate o-spanins(s): "+str(len(v))+"\n")
                    if each_spanin_type == "SEPARATED":
                        for k, v in SPANIN_TYPES['SEPARATED'].items():
                            f.writelines(""+str(k)+ " ==> Amount of corresponding candidate o-spanins(s): "+str(len(v))+"\n")
                    if each_spanin_type == "OVERLAPPED":
                        for k, v in SPANIN_TYPES['OVERLAPPED'].items():
                            f.writelines(""+str(k)+" ==> Amount of corresponding candidate o-spanins(s): "+str(len(v))+"\n")
            except TypeError:
                continue
        f.write("\n=~~~~~= Tally from ALL spanin types =~~~~~=\n")
        f.writelines("Total Candidates = "+str(spanins['total_amount'])+"\n")
        f.writelines("Total Unique Candidate Pairs = "+str(spanins['total_unique'])+"\n")

    args.putative_isp_fasta_file = open(args.putative_isp_fasta_file.name, "r")
    isp_full = tuple_fasta(args.putative_isp_fasta_file)

    args.putative_osp_fasta_file = open(args.putative_osp_fasta_file.name, "r")
    osp_full = tuple_fasta(args.putative_osp_fasta_file)

    #print(isp_full)
    isp_seqs = []
    osp_seqs = []
    for isp_tupe in isp_full:
        #print(isp_tupe)
        for pisp, posp in embedded.items():
            #print(f"ISP = searching for {pisp} in {isp_tupe[0]}")
            if re.search(("("+str(pisp)+")\D"), isp_tupe[0]):
                #print(isp_tupe[0])
                #print(peri_count)
                peri_count = str.split(isp_tupe[0],"~=")[1]
                isp_seqs.append((pisp,isp_tupe[1],peri_count))
    #print(isp_seqs)
    for osp_tupe in osp_full:
        for pisp, posp in embedded.items():
            for data in posp:
                #print(f"OSP = searching for {data[3]} in {osp_tupe[0]}, coming from this object: {data}")
                if re.search(("("+str(data[3])+")\D"), osp_tupe[0]):
                    peri_count = str.split(osp_tupe[0],"~=")[1]
                    osp_seqs.append((data[3],osp_tupe[1],peri_count))

    with args.embedded_txt as f:
        f.write("================ embedded spanin candidates =================\n")
        f.write("isp\tisp_start\tisp_end\tosp\tosp_start\tosp_end\tstrand\tpair_number\n")
        if embedded != {}:
            #print(embedded)
            for pisp, posp in embedded.items():
                #print(f"{pisp} - {posp}")
                f.write(pisp + "\n")
                for each_posp in posp:
                    #print(posp)
                    f.write(
                        "\t{}\t{}\t{}\t{}\t{}\t{}\t".format(
                            each_posp[1],
                            each_posp[2],
                            each_posp[3],
                            each_posp[4],
                            each_posp[5],
                            each_posp[6],
                        )
                    )
                    if each_posp[6] == "+":
                        if each_posp[2] in pair_dict['pairs']['pair_number'].keys():
                            f.write(""+str(pair_dict['pairs']['pair_number'][each_posp[2]])+"\n")
                    elif each_posp[6] == "-":
                        if each_posp[1] in pair_dict['pairs']['pair_number'].keys():
                            f.write(""+str(pair_dict['pairs']['pair_number'][each_posp[1]])+"\n")
        else:
            f.write("nothing found")

    with open(args.embedded_txt.name, "a") as f:
        f.write("\n================= embedded candidate sequences ================\n")
        f.write("======================= isp ==========================\n\n")
        for isp_data in isp_seqs:
            #print(isp_data)
            f.write(">isp_orf::{}-peri_count~={}\n{}\n".format(isp_data[0],isp_data[2],lineWrapper(isp_data[1])))
        f.write("\n======================= osp ========================\n\n")
        for osp_data in osp_seqs:
            f.write(">osp_orf::{}-peri_count~={}\n{}\n".format(osp_data[0],osp_data[2],lineWrapper(osp_data[1])))

    args.putative_isp_fasta_file = open(args.putative_isp_fasta_file.name, "r")
    isp_full = tuple_fasta(args.putative_isp_fasta_file)

    args.putative_osp_fasta_file = open(args.putative_osp_fasta_file.name, "r")
    osp_full = tuple_fasta(args.putative_osp_fasta_file)

    isp_seqs = []
    osp_seqs = []
    for isp_tupe in isp_full:
        peri_count = str.split(isp_tupe[0],"~=")[1]
        for pisp, posp in overlap.items():
            if re.search(("("+str(pisp)+")\D"), isp_tupe[0]):
                peri_count = str.split(isp_tupe[0],"~=")[1]
                isp_seqs.append((pisp,isp_tupe[1],peri_count))

    for osp_tupe in osp_full:
        for pisp, posp in overlap.items():
            for data in posp:
                if re.search(("("+str(data[3])+")\D"), osp_tupe[0]):
                    peri_count = str.split(osp_tupe[0],"~=")[1]
                    osp_seqs.append((data[3],osp_tupe[1],peri_count))


    
    with args.overlap_txt as f:
        f.write("================ overlap spanin candidates =================\n")
        f.write("isp\tisp_start\tisp_end\tosp\tosp_start\tosp_end\tstrand\tpair_number\n")
        if overlap != {}:
            for pisp, posp in overlap.items():
                f.write(pisp + "\n")
                for each_posp in posp:
                    f.write(
                        "\t{}\t{}\t{}\t{}\t{}\t{}\t".format(
                            each_posp[1],
                            each_posp[2],
                            each_posp[3],
                            each_posp[4],
                            each_posp[5],
                            each_posp[6],
                        )
                    )
                    if each_posp[6] == "+":
                        if each_posp[2] in pair_dict['pairs']['pair_number'].keys():
                            #print('ovl ; +')
                            f.write(""+str(pair_dict['pairs']['pair_number'][each_posp[2]])+"\n")
                    elif each_posp[6] == "-":
                        if each_posp[1] in pair_dict['pairs']['pair_number'].keys():
                            f.write(""+str(pair_dict['pairs']['pair_number'][each_posp[1]])+"\n")
        else:
            f.write("nothing found")

    with open(args.overlap_txt.name, "a") as f:
        #print(isp_seqs)
        f.write("\n================= overlap candidate sequences ================\n")
        f.write("======================= isp ==========================\n\n")
        for isp_data in isp_seqs:
            f.write(">isp_orf::{}-pericount~={}\n{}\n".format(isp_data[0],isp_data[2],lineWrapper(isp_data[1])))
        f.write("\n======================= osp ========================\n\n")
        for osp_data in osp_seqs:
            f.write(">osp_orf::{}-pericount~={}\n{}\n".format(osp_data[0],osp_data[2],lineWrapper(osp_data[1])))

    args.putative_isp_fasta_file = open(args.putative_isp_fasta_file.name, "r")
    isp_full = tuple_fasta(args.putative_isp_fasta_file)
    args.putative_osp_fasta_file = open(args.putative_osp_fasta_file.name, "r")
    osp_full = tuple_fasta(args.putative_osp_fasta_file)

    isp_seqs = []
    osp_seqs = []
    for isp_tupe in isp_full:
        for pisp, posp in separate.items():
            if re.search(("("+str(pisp)+")\D"), isp_tupe[0]):
                peri_count = str.split(isp_tupe[0],"~=")[1]
                isp_seqs.append((pisp,isp_tupe[1],peri_count))
    #print(isp_seqs)
    for osp_tupe in osp_full:
        for pisp, posp in separate.items():
            for data in posp:
                if re.search(("("+str(data[3])+")\D"), osp_tupe[0]):
                    peri_count = str.split(osp_tupe[0],"~=")[1]
                    osp_seqs.append((data[3],osp_tupe[1],peri_count))

    with args.separate_txt as f:
        f.write("================ separated spanin candidates =================\n")
        f.write("isp\tisp_start\tisp_end\tosp\tosp_start\tosp_end\tstrand\tpair_number\n")
        if separate != {}:
            for pisp, posp in separate.items():
                f.write(pisp + "\n")
                for each_posp in posp:
                    f.write(
                        "\t{}\t{}\t{}\t{}\t{}\t{}\t".format(
                            each_posp[1],
                            each_posp[2],
                            each_posp[3],
                            each_posp[4],
                            each_posp[5],
                            each_posp[6],
                        )
                    )
                    if each_posp[6] == "+":
                        if each_posp[2] in pair_dict['pairs']['pair_number'].keys():
                            f.write(""+str(pair_dict['pairs']['pair_number'][each_posp[2]])+"\n")
                    elif each_posp[6] == "-":
                        if each_posp[1] in pair_dict['pairs']['pair_number'].keys():
                            f.write(""+str(pair_dict['pairs']['pair_number'][each_posp[1]])+"\n")
        else:
            f.write("nothing found")

    with open(args.separate_txt.name, "a") as f:
        f.write("\n================= separated candidate sequences ================\n")
        f.write("======================= isp ==========================\n\n")
        for isp_data in isp_seqs:
            f.write(">isp_orf::{}-pericount~={}\n{}\n".format(isp_data[0],isp_data[2],lineWrapper(isp_data[1])))
        f.write("\n======================= osp ========================\n\n")
        for osp_data in osp_seqs:
            f.write(">osp_orf::{}-pericount~={}\n{}\n".format(osp_data[0],osp_data[2],lineWrapper(osp_data[1])))
