
def fasta_from_SAR_dict(sar_dict,fa_file):
    """ makes a multi fasta with candidates from SAR dictionary """
    with fa_file as f:
        for data in sar_dict.values():
            f.writelines(">{}\n".format(data["description"]))
            f.writelines("{}\n".format(data["sequence"]))

def stat_file_from_SAR_dict(sar_dict, stat_file, tmd_min, tmd_max):
    """ summary statistics from SAR finder function """
    with stat_file as f:
        f.writelines("..........:::::: Candidate SAR Proteins ::::::..........\n\n")
        if sar_dict:
            for data in sar_dict.values():
                f.writelines("Protein Description and Name: {}\n".format(data["description"]))
                f.writelines("Protein Sequence: {}\n".format(data["sequence"]))
                f.writelines("Protein Length: {}\n".format(data["size"]))
                f.writelines("Percent G and A content: {}\n".format(data["GAcont"]))
                f.writelines("SAR Criteria matching region(s)\n")
                for tmd_size in range(tmd_min, tmd_max, 1):
                    if "TMD_"+str(tmd_size) in data:
                        f.writelines("TMD length of {}:\n".format(tmd_size))
                        for each_match in data["TMD_"+str(tmd_size)]:
                            f.writelines("TMD domain sequence: {}\n".format(each_match[0]))
                            f.writelines("N-term sequence: {}\n".format(each_match[1]))
                            f.writelines("N-term net charge: {}\n".format(each_match[2]))
                f.writelines("========================================================\n\n")
        else:
            f.writelines("No candidate SAR Proteins found")


if __name__ == "__main__":
    pass