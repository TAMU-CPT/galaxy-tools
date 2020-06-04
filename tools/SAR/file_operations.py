
def fasta_from_SAR_dict(sar_dict,fa_file):
    """ makes a multi fasta with candidates from SAR dictionary """
    with fa_file as f:
        for data in sar_dict.values():
            f.writelines(">{}\n".format(data["description"]))
            f.writelines("{}\n".format(data["sequence"]))

def gff3_from_SAR_dict(sar_dict,gff3_file):
    """ make a multi gff3 with candidates from SAR dictionary """
    pass

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
                            f.writelines("\nPotential SAR domain sequence: {}\n".format(each_match[0]))
                            f.writelines("N-term sequence: {}\n".format(each_match[1]))
                            f.writelines("N-term net charge: {}\n".format(each_match[2]))
                            for each_perc_calc in each_match[3]:
                                f.writelines("Percent {} content: {}%\n".format(each_perc_calc[0],each_perc_calc[1]))
                            f.writelines("N-term coords: {}\n".format(each_match[4]))
                            f.writelines("SAR coords: {}\n".format(each_match[5]))
                            f.writelines("C-term coords: {}\n".format(each_match[6]))
                            f.writelines("SAR start: {}\n".format(each_match[7]))
                f.writelines("========================================================\n\n")
        else:
            f.writelines("No candidate SAR Proteins found")


if __name__ == "__main__":
    pass