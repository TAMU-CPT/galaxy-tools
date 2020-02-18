##### findSpanin.pl --> findSpanin.py
######### Incooperated from the findSpanin.pl script, but better and more snakey.

import argparse
import os
from spaninFuncs import getDescriptions, grabLocs, spaninProximity, splitStrands

### Requirement Inputs
#### INPUT : Genomic FASTA
#### PARAMETERS :

###############################################################################

if __name__ == '__main__':

    # Common parameters for both ISP / OSP portion of script

    parser = argparse.ArgumentParser(description='Trim the putative protein candidates and find potential i-spanin / o-spanin pairs')

    parser.add_argument('putative_isp_fasta_file', type=argparse.FileType("r"), 
                        help='Putative i-spanin Fasta file, output of "generate-putative-isp"') # the "input" argument

    parser.add_argument('putative_osp_fasta_file', type=argparse.FileType("r"), 
                        help='Putative o-spanin Fasta file, output of "generate-putative-osp"')

    parser.add_argument('--max_isp_osp_distance', dest='max_isp_osp_distance', default=10, help='max distance from end of i-spanin to start of o-spanin, measured in AAs')

    parser.add_argument('--strand', dest='strand', default='+', help='strand to investigate matches, + or -')
    parser.add_argument('-v', action='version', version='0.3.0') # Is this manually updated?
    args = parser.parse_args()

    
    isp = getDescriptions(args.putative_isp_fasta_file.name)
    print(len(isp))
    osp = getDescriptions(args.putative_osp_fasta_file.name)
    
    strand_isp = []
    strand_osp = []
    for desc in isp: # will retrieve only + or - strand for analysis
        text = splitStrands(desc,args.strand)
        strand_isp.append(text)
    for desc in osp:
        text = splitStrands(desc,args.strand)
        strand_osp.append(text)

    strand_isp = [i for i in strand_isp if i] # filtering out Nones
    strand_osp = [ii for ii in strand_osp if ii] # filtering out Nones

    data_isp = []
    data_osp = []
    for desc in strand_isp:
        d = grabLocs(desc)
        data_isp.append(d)
    
    for desc in strand_osp:
        d = grabLocs(desc)
        data_osp.append(d)
    
    print(data_isp)
    print(data_osp)

    embedded, overlap, upstream = spaninProximity(data_isp,
                                                    data_osp, 
                                                    max_dist=args.max_isp_osp_distance*3)
    print('------candidates-------')
    #print('\n')
    print(len(embedded))
    #print('\n')
    print(len(overlap))
    #print('\n')
    print(len(upstream))







