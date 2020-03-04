#!/usr/bin/env python
import numpy
import argparse
from Bio import SeqIO
import logging

logging.basicConfig(level=logging.INFO)


class ABIFParser(object):
    def __init__(self, abif_file):
        records = list(SeqIO.parse(abif_file, "abi"))

        for record in records:
            self.abif = record.annotations["abif_raw"]

        # Copy+pasted from seqtrace-0.9.0
        #
        # According to the ABIF documentation, ABIF files (after base calling)
        # should contain two base call entries ("PBAS"): one containing
        # "sequence characters edited by user" (entry number 1), and one
        # containing "sequence characters as called by Basecaller" (entry
        # number 2). These two entries will, in most cases, contain identical
        # sequence data. This method follows the same convention used by the
        # Staden package (see seqIOABI.c), which is to only look at entry 1
        # (the user-edited sequence) and ignore entry 2.
        self.basecalls = self.abif["PBAS1"]

        # There is an inconsistency in the ABIF file format documentation
        # regarding the data format of the confidence scores. The data format
        # ID (as actually found in a .ab1 file) is 2, indicating the values are
        # signed 1-byte integers. The ABIF documentation, however, sugggests
        # the values can range from 0-255 (i.e., an unsigned 1-byte integer).
        # In practice, the actual values do not appear to exceed 61, making the
        # distinction between signed/unsigned irrelevant. For now, the data
        # format ID is taken as the correct indication of the underlying data
        # format.
        #
        # According to the ABIF documentation, ABIF files (after base calling)
        # should contain two quality value (QV) entries ("PCON"): one
        # containing QVs "as edited by user" (entry number 1), and one
        # containing QVs "as called by Basecaller" (entry number 2). These two
        # entries will, in most cases, contain identical values. This method
        # follows the same convention used by the Staden package (see
        # seqIOABI.c), which is to only look at entry 1 (the user-edited QVs)
        # and ignore entry 2.
        self.bcconf = map(ord, self.abif["PCON1"])

        self.base_order = self.abif["FWO_1"]  # filter wheel order

        # This is the ID for the first 'DATA' index entry that points to the
        # processed trace data. The man page for the Staden program
        # convert_trace states that IDs 9-12 contain the processed data; IDs
        # 1-4 contain the raw data. The ABIF documentation from ABI also
        # suggests that IDs 1-8 will always contain raw data, and 9-12 will
        # contain the processed data. Is this always correct?

        sn = {}
        for i, char in enumerate(list(self.abif["FWO_1"])):
            sn[char] = {
                "wavelength": self.abif["DyeW" + str(i + 1)],
                "correction": self.abif["S/N%1"][i],
                "data": numpy.array(self.abif["DATA" + str(i + 9)]).astype(numpy.float),
            }
        print(list(map(ord, self.abif["PCON2"])))
        print(len(sn["A"]["data"]))

        spacing = self.abif["SPAC1"]
        if spacing < 0:
            spacing = float(self.basepos[-1] - self.basepos[0]) / (
                len(self.basepos) - 1
            )
        print(spacing)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot AB1 file")
    parser.add_argument("abif_file", type=argparse.FileType("rb"), help="ABIF/AB1 file")
    args = parser.parse_args()
    ap = ABIFParser(args.abif_file)
    # plot(**vars(args))
