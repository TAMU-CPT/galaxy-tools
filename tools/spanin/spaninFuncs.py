"""
PREMISE
### Functions/Classes that are used in both generate-putative-osp.py and generate-putative-isp.py
###### Main premise here is to make the above scripts a little more DRY, as well as easily readable for execution.
###### Documentation will ATTEMPT to be thourough here
"""

import re
from Bio import SeqIO
from Bio import Seq

# Not written in OOP for a LITTLE bit of trying to keep the complication down in case adjustments are needed by someone else.
# Much of the manipulation is string based; so it should be straightforward as well as moderately quick
################## GLobal Variables
Lys = "K"


def check_back_end_snorkels(seq, tmsize):
    """
        Searches through the backend of a potential TMD snorkel. This is the 2nd part of a TMD snorkel lysine match.
        --> seq : should be the sequence fed from the "search_region" portion of the sequence
        --> tmsize : size of the potential TMD being investigated
    """
    found = []
    if seq[tmsize - 4] == Lys and re.search(("[FIWLVMYCATGS]"), seq[tmsize - 5]):
        found = "match"
        return found
    elif seq[tmsize - 3] == Lys and re.search(("[FIWLVMYCATGS]"), seq[tmsize - 4]):
        found = "match"
        return found
    elif seq[tmsize - 2] == Lys and re.search(("[FIWLVMYCATGS]"), seq[tmsize - 3]):
        found = "match"
        return found
    elif seq[tmsize - 1] == Lys and re.search(("[FIWLVMYCATGS]"), seq[tmsize - 2]):
        found = "match"
        return found
    else:
        found = "NOTmatch"
        return found


def prep_a_gff3(fa, spanin_type):
    """
        Function parses an input detailed 'fa' file and outputs a 'gff3' file
        ---> fa = input .fa file
        ---> output = output a returned list of data, easily portable to a gff3 next
        ---> spanin_type = 'isp' or 'osp'
    """
    fa_zip = tuple_fasta(fa)
    data = []
    for a_pair in fa_zip:
        # print(a_pair)
        if re.search(("(\[1\])"), a_pair[0]):
            strand = "+"
        elif re.search(("(\[-1\])"), a_pair[0]):
            strand = "-"  # column 7
        start = re.search(("[\d]+\.\."), a_pair[0]).group(0).split("..")[0]  # column 4
        end = re.search(("\.\.[\d]+"), a_pair[0]).group(0).split("..")[1]  # column 5
        orfid = re.search(("(ORF)[\d]+"), a_pair[0]).group(0)  # column 1
        if spanin_type == "isp":
            methodtype = "CDS"  # column 3
            spanin = "isp"
        elif spanin_type == "osp":
            methodtype = "CDS"  # column 3
            spanin = "osp"
        else:
            print("need to input spanin type")
            break
        source = "cpt.py|putative-*.py"  # column 2
        score = "."  # column 6
        phase = "."  # column 8
        seq = a_pair[1] + ";Alias=" + spanin  # column 9
        sequence = [[orfid, source, methodtype, start, end, score, strand, phase, seq]]
        data += sequence
    return data


def write_gff3(data, output="results.gff3"):
    """
        Parses results from prep_a_gff3 into a gff3 file
        ---> input : list from prep_a_gff3
        ---> output : gff3 file
    """
    data = data
    filename = output
    with filename as f:
        f.write("#gff-version 3\n")
        for value in data:
            f.write(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    value[0],
                    value[1],
                    value[2],
                    value[3],
                    value[4],
                    value[5],
                    value[6],
                    value[7],
                    value[8],
                )
            )
    f.close()


def find_tmd(pair, minimum=10, maximum=30, TMDmin=10, TMDmax=20):
    """ 
        Function that searches for lysine snorkels and then for a spanning hydrophobic region that indicates a potential TMD
        ---> pair : Input of tuple with description and AA sequence (str)
        ---> minimum : How close from the initial start codon a TMD can be within
        ---> maximum : How far from the initial start codon a TMD can be within
        ---> TMDmin : The minimum size that a transmembrane can be (default = 10)
        ---> TMDmax : The maximum size tha ta transmembrane can be (default = 20)
    """
    # hydrophobicAAs = ['F', 'I', 'W', 'L', 'V', 'M', 'Y', 'C', 'A', 'T', 'G', 'S']
    tmd = []
    s = str(pair[1])  # sequence being analyzed
    # print(s) # for trouble shooting
    search_region = s[minimum - 1 : maximum + 1]
    # print(search_region) # for trouble shooting

    for tmsize in range(TMDmin, TMDmax + 1, 1):
        # print('==============='+str(tmsize)+'================') # print for troubleshooting
        pattern = (
            "['FIWLVMYCATGS']{" + str(tmsize) + "}"
        )  # searches for these hydrophobic residues tmsize total times
        if re.search(
            ("[K]"), search_region[1:8]
        ):  # grabbing one below with search region, so I want to grab one ahead here when I query.
            store_search = re.search(
                ("[K]"), search_region[1:8]
            )  # storing regex object
            where_we_are = store_search.start()  # finding where we got the hit
            if re.search(
                ("[FIWLVMYCATGS]"), search_region[where_we_are + 1]
            ) and re.search(
                ("[FIWLVMYCATGS]"), search_region[where_we_are - 1]
            ):  # hydrophobic neighbor
                try:
                    backend = check_back_end_snorkels(search_region, tmsize)
                    if backend == "match":
                        tmd.append(pair)
                    else:
                        continue
                except (IndexError, TypeError):
                    continue
            else:
                continue
        elif re.search((pattern), search_region):
            try:
                tmd.append(pair)
            except (IndexError, TypeError):
                continue
        else:
            continue

        return tmd


def find_lipobox(pair, minimum=10, maximum=30, min_after=10, max_after=50, regex=1):
    """
        Function that takes an input tuple, and will return pairs of sequences to their description that have a lipoobox
        ---> minimum - min distance from start codon to first AA of lipobox
        ---> maximum - max distance from start codon to first AA of lipobox
        ---> regex - option 1 (default) => more strict regular expression ; option 2 => looser selection, imported from LipoRy
        
    """
    if regex == 1:
        pattern = "[ILMFTV][^REKD][GAS]C"  # regex for Lipobox from findSpanin.pl
    elif regex == 2:
        pattern = "[ACGSILMFTV][^REKD][GAS]C"  # regex for Lipobox from LipoRy

    candidates = []
    s = str(pair[1])
    # print(s) # trouble shooting
    search_region = s[minimum-1 : maximum + 5] # properly slice the input... add 4 to catch if it hangs off at max input
    # print(search_region) # trouble shooting
    # for each_pair in pair:
    # print(s)
    if re.search((pattern), search_region):  # lipobox must be WITHIN the range...
        # searches the sequence with the input RegEx AND omits if
        #g = re.search((pattern), search_region).group() # find the exact group match
        #if min_after < len(s) - re.search((g), s).end() < max_after: # find the lipobox end region
        candidates.append(pair)
        # print('passed') # trouble shooting
        return candidates
    else:
        # print('didnotpass') # trouble shooting
        pass


def tuple_fasta(fasta_file):
    """
        #### INPUT: Fasta File
        #### OUTPUT: zipped (zip) : pairwise relationship of description to sequence
        #### 
    """
    fasta = SeqIO.parse(fasta_file, "fasta")
    descriptions = []
    sequences = []
    for r in fasta:  # iterates and stores each description and sequence
        description = r.description
        sequence = str(r.seq)
        if (
            sequence[0] != "I"
        ):  # the translation table currently has I as a potential start codon ==> this will remove all ORFs that start with I
            descriptions.append(description)
            sequences.append(sequence)
        else:
            continue

    return zip(descriptions, sequences)


def lineWrapper(text, charactersize=60):

    if len(text) <= charactersize:
        return text
    else:
        return (
            text[:charactersize]
            + "\n"
            + lineWrapper(text[charactersize:], charactersize)
        )


def getDescriptions(fasta):
    """
    Takes an output FASTA file, and parses retrieves the description headers. These headers contain information needed
    for finding locations of a potential i-spanin and o-spanin proximity to one another.
    """
    desc = []
    with fasta as f:
        for line in f:
            if line.startswith(">"):
                desc.append(line)
    return desc


def splitStrands(text, strand="+"):
    # positive_strands = []
    # negative_strands = []
    if strand == "+":
        if re.search(("(\[1\])"), text):
            return text
    elif strand == "-":
        if re.search(("(\[-1\])"), text):
            return text
    # return positive_strands, negative_strands


def parse_a_range(pair, start, end):
    """
        Takes an input data tuple from a fasta tuple pair and keeps only those within the input sequence range
        ---> data : fasta tuple data
        ---> start : start range to keep
        ---> end : end range to keep (will need to + 1)
    """
    matches = []
    for each_pair in pair:

        s = re.search(("[\d]+\.\."), each_pair[0]).group(0)  # Start of the sequence
        s = int(s.split("..")[0])
        e = re.search(("\.\.[\d]+"), each_pair[0]).group(0)
        e = int(e.split("..")[1])
        if start - 1 <= s and e <= end + 1:
            matches.append(each_pair)
        else:
            continue
    # else:
    # continue
    # if matches != []:
    return matches
    # else:
    # print('no candidates within selected range')


def grabLocs(text):
    """
    Grabs the locations of the spanin based on NT location (seen from ORF). Grabs the ORF name, as per named from the ORF class/module
    from cpt.py
    """
    start = re.search(("[\d]+\.\."), text).group(
        0
    )  # Start of the sequence ; looks for [numbers]..
    end = re.search(("\.\.[\d]+"), text).group(
        0
    )  # End of the sequence ; Looks for ..[numbers]
    orf = re.search(("(ORF)[\d]+"), text).group(
        0
    )  # Looks for ORF and the numbers that are after it

    start = int(start.split("..")[0])
    end = int(end.split("..")[1])

    vals = [start, end, orf]

    """
    store_vals = []
    for r in vals:
        if r is not None:
            store_vals.append(r.group(0))
    """
    return vals


def spaninProximity(isp, osp, max_dist=30, strand="+"):
    """
    _NOTE THIS FUNCTION COULD BE MODIFIED TO RETURN SEQUENCES_
    Compares the locations of i-spanins and o-spanins. max_dist is the distance in NT measurement from i-spanin END site
    to o-spanin START. The user will be inputting AA distance, so a conversion will be necessary (<user_input> * 3)
    INPUT: list of OSP and ISP candidates
    OUTPUT: Return (improved) candidates for overlapping, embedded, and separate list
    """
    if strand == "+":
        embedded = {}
        overlap = {}
        separate = {}
        for iseq in isp:
            embedded[iseq[2]] = []
            overlap[iseq[2]] = []
            separate[iseq[2]] = []
            # print(iseq)
            for oseq in osp:
                # print(oseq)
                if iseq[0] < oseq[0] < iseq[1] and oseq[1] < iseq[1]:
                    ### EMBEDDED ###
                    combo = [
                        iseq[0],
                        iseq[1],
                        oseq[2],
                        oseq[0],
                        oseq[1],
                    ]  # ordering a return for dic
                    embedded[iseq[2]] += [combo]
                elif iseq[0] < oseq[0] <= iseq[1] and oseq[1] > iseq[1]:
                    ### OVERLAP / SEPARATE ###
                    if (iseq[1] - oseq[0]) < 6:
                        combo = [iseq[0], iseq[1], oseq[2], oseq[0], oseq[1]]
                        separate[iseq[2]] += [combo]
                    else:
                        combo = [iseq[0], iseq[1], oseq[2], oseq[0], oseq[1]]
                        overlap[iseq[2]] += [combo]
                elif iseq[1] <= oseq[0] <= iseq[1] + max_dist:
                    combo = [iseq[0], iseq[1], oseq[2], oseq[0], oseq[1]]
                    separate[iseq[2]] += [combo]
                else:
                    continue
    elif strand == "-":
        embedded = {}
        overlap = {}
        separate = {}
        for iseq in isp:
            embedded[iseq[2]] = []
            overlap[iseq[2]] = []
            separate[iseq[2]] = []
            for oseq in osp:
                if iseq[0] <= oseq[1] <= iseq[1] and oseq[0] > iseq[0]:
                    ### EMBEDDED ###
                    combo = [
                        iseq[0],
                        iseq[1],
                        oseq[2],
                        oseq[0],
                        oseq[1],
                    ]  # ordering a return for dict
                    embedded[iseq[2]] += [combo]
                elif iseq[0] <= oseq[1] <= iseq[1] and oseq[0] < iseq[0]:
                    if (oseq[1] - iseq[0]) < 6:
                        combo = [iseq[0], iseq[1], oseq[2], oseq[0], oseq[1]]
                        separate[iseq[2]] += [combo]
                    else:
                        combo = [iseq[0], iseq[1], oseq[2], oseq[0], oseq[1]]
                        overlap[iseq[2]] += [combo]
                elif iseq[0] - 10 < oseq[1] < iseq[0]:
                    combo = [iseq[0], iseq[1], oseq[2], oseq[0], oseq[1]]
                    separate[iseq[2]] += [combo]
                else:
                    continue
    else:
        print("please insert a strand")
        pass

    embedded = {k: embedded[k] for k in embedded if embedded[k]}
    overlap = {k: overlap[k] for k in overlap if overlap[k]}
    separate = {k: separate[k] for k in separate if separate[k]}
    return embedded, overlap, separate


############################################### TEST RANGE #########################################################################
####################################################################################################################################
if __name__ == "__main__":

    #### TMD TEST
    test_desc = ["one", "two", "three", "four", "five"]
    test_seq = [
        "XXXXXXXXXXXXXXXFMCFMCFMCFMCFMCXXXXXXXXXXXXXXXXXXXXXXXXXX",
        "XXXXXXXXAAKKKKKKKKKKKKKKKXXXXXXXXXXXXX",
        "XXXXXXX",
        "XXXXXXXXXXXKXXXXXXXXXX",
        "XXXXXXXXXXAKXXXXXXXXXXAKXXXXXXXX",
    ]
    # for l in
    # combo = zip(test_desc,test_seq)
    pairs = zip(test_desc, test_seq)
    tmd = []
    for each_pair in pairs:
        # print(each_pair)
        try:
            tmd += find_tmd(pair=each_pair)
        except (IndexError, TypeError):
            continue
        # try:s = each_pair[1]
        # tmd += find_tmd(seq=s, tmsize=15)
    # print('\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    # print(tmd)
    # print('\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')

    #### tuple-fasta TEST
    # fasta_file = 'out_isp.fa'
    # ret = tuple_fasta(fasta_file)
    # print('=============')
    # for i in ret:
    # print(i[1])

    #### LipoBox TEST
    test_desc = ["one", "two", "three", "four", "five", "six", "seven"]
    test_seq = [
        "XXXXXXXXXTGGCXXXXXXXXXXXXXXXX",
        "XXXXXXXXAAKKKKKKKKKKKKKKKXXXXXXXXXXXXX",
        "XXXXXXX",
        "AGGCXXXXXXXXXXXXXXXXXXXXTT",
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXTGGC",
        "XXXXXXXXXXXXXXXXXXXXXXXXXXTGGC",
        "MSTLRELRLRRALKEQSMRYLLSIKKTLPRWKGALIGLFLICVATISGCASESKLPEPPMVSVDSSLMVEPNLTTEMLNVFSQ*",
    ]
    pairs = zip(test_desc, test_seq)
    lipo = []
    for each_pair in pairs:
        print(each_pair)
        # try:
        try:
            lipo += find_lipobox(pair=each_pair, regex=2)  # , minimum=8)
        except TypeError:  # catches if something doesnt have the min/max requirements (something is too small)
            continue
        # except:
        # continue
    # print('\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    #############################3
    # g = prep_a_gff3(fa='putative_isp.fa', spanin_type='isp')
    # print(g)
    # write_gff3(data=g)
