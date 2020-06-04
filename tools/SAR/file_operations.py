
import pandas as pd


def fasta_from_SAR_dict(sar_dict,fa_file):
    """ makes a multi fasta with candidates from SAR dictionary """
    with fa_file as f:
        for data in sar_dict.values():
            f.writelines(">{}\n".format(data["description"]))
            f.writelines("{}\n".format(data["sequence"]))

def gff3_from_SAR_dict(sar_dict,gff3_file):
    """ make a multi gff3 with candidates from SAR dictionary """
    gff3_cols = ["Seqid","Source","Type","Start","End","Score","Strand","Phase","Attributes"]
    with gff3_file as f:
        f.writelines(f"{gff3_cols[0]}\t{gff3_cols[1]}\t{gff3_cols[2]}\t{gff3_cols[3]}\t{gff3_cols[4]}\t{gff3_cols[5]}\t{gff3_cols[6]}\t{gff3_cols[7]}\t{gff3_cols[8]}\n")
        if sar_dict:
            #print(sar_dict)
            for name, data in sar_dict.items():
                if len(data["TMD_"+str(data["biggest_sar"])]) > 1: # might need to be ["size"][-1]...dont remember if this will be ordered in reverse or not...
                    values = []
                    for idx, value in enumerate(data["TMD_"+str(data["biggest_sar"])][0]):
                        print(value)
                        values.append([idx,value])
                    print(values)
                    min_idx = min(values[1])
                else:
                    min_idx = 0
                f.writelines("##gff-version 3\n")
                f.writelines(f"##sequence-region {name}\n")
                n_start, n_end = split_seq_string(data["TMD_"+str(data["biggest_sar"])][min_idx][4])
                sar_start, sar_end = split_seq_string(data["TMD_"+str(data["biggest_sar"])][min_idx][5])
                c_start, c_end = split_seq_string(data["TMD_"+str(data["biggest_sar"])][min_idx][6])
                f.writelines(f'{name}\tSAR_finder\tTopological domain\t{n_start}\t{n_end}\t.\t.\t.\tNote=N-terminus Charge is {data["TMD_"+str(data["biggest_sar"])][min_idx][2]}\n')
                f.writelines(f'{name}\tSAR_finder\tSAR domain\t{sar_start}\t{sar_end}\t.\t.\t.\tNote=%{[perc for perc in data["TMD_"+str(data["biggest_sar"])][min_idx][3]]}\n')
                f.writelines(f'{name}\tSAR_finder\tTopological domain\t{c_start}\t{c_end}\t.\t.\t.\tNote=C-terminus\n')
        else:
            f.writelines("##gff-version 3\n")
            f.writelines(f"##sequence-region {name}\n")


def tab_from_SAR_dict(sar_dict,stat_file,hydrophillic_res, sar_min, sar_max):
    """ convert SAR dict to a dataframe """
    columns = ["Name","Protein Sequence","Protein Length","SAR Length","Putative SAR Sequence","SAR Start Location",[f"{res}%" for res in hydrophillic_res],"N-term Sequence","N-term net Charge"]
    with stat_file as f:
        f.writelines(f"{columns[0]}\t{columns[1]}\t{columns[2]}\t{columns[3]}\t{columns[4]}\t{columns[5]}\t{columns[6]}\t{columns[7]}\t{columns[8]}\n")
        if sar_dict:
            for name, data in sar_dict.items():
                for tmd_size in range(sar_max, sar_min-1, -1):
                    if "TMD_"+str(tmd_size) in data:
                        for each_match in data["TMD_"+str(tmd_size)]:
                            if each_match != [""]:
                                f.writelines(f'{name}\t{data["sequence"]}\t{data["size"]}\t{tmd_size}\t{each_match[0]}\t{int(each_match[-1])+1}\t{[perc[1] for perc in each_match[3]]}\t{each_match[1]}\t{each_match[2]}\n')
                            else:
                                continue

def stat_file_from_SAR_dict(sar_dict, stat_file, sar_min, sar_max):
    """ summary statistics from SAR finder function """
    with stat_file as f:
        f.writelines("..........:::::: Candidate SAR Proteins ::::::..........\n\n")
        if sar_dict:
            for data in sar_dict.values():
                f.writelines("Protein Description and Name: {}\n".format(data["description"]))
                f.writelines("Protein Sequence: {}\n".format(data["sequence"]))
                f.writelines("Protein Length: {}\n".format(data["size"]))
                f.writelines("SAR Criteria matching region(s)\n")
                for tmd_size in range(sar_max, sar_min-1, -1):
                    if "TMD_"+str(tmd_size) in data:
                        f.writelines("\nSAR length of {}:\n".format(tmd_size))
                        for each_match in data["TMD_"+str(tmd_size)]:
                            if each_match != ['']:
                                f.writelines("\nPotential SAR domain sequence: {}\n".format(each_match[0]))
                                f.writelines("N-term sequence: {}\n".format(each_match[1]))
                                f.writelines("N-term net charge: {}\n".format(each_match[2]))
                                for each_perc_calc in each_match[3]:
                                    f.writelines("Percent {} content: {}%\n".format(each_perc_calc[0],each_perc_calc[1]))
                                f.writelines("N-term coords: {}\n".format(each_match[4]))
                                f.writelines("SAR coords: {}\n".format(each_match[5]))
                                f.writelines("C-term coords: {}\n".format(each_match[6]))
                                f.writelines("SAR start: {}\n".format(each_match[7]))
                            else:
                                continue
                f.writelines("========================================================\n\n")
        else:
            f.writelines("No candidate SAR Proteins found")

def split_seq_string(input_range, python_indexing=True):
    """ splits a #..# sequence into the two respective starts and ends, if python indexing, adds 1, otherwise keeps """
    if python_indexing:
        values = input_range.split("..")
        start =int(values[0]) + 1
        end = int(values[1]) + 1
    else:
        values = input_range.split("..")
        start = values[0]
        end = values[1]

    return start, end

if __name__ == "__main__":
    pass