##### findSpanin.pl --> findSpanin.py
######### Much of this code is very "blocked", in the sense that one thing happens...then a function happens on the return...then another function...etc...etc...

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
                        help='Putative i-spanin FASTA file, output of "generate-putative-isp"') # the "input" argument

    parser.add_argument('putative_osp_fasta_file', type=argparse.FileType("r"), 
                        help='Putative o-spanin FASTA file, output of "generate-putative-osp"')

    parser.add_argument('--max_isp_osp_distance', dest='max_isp_osp_distance', default=10, type=int, help='max distance from end of i-spanin to start of o-spanin, measured in AAs')

    parser.add_argument('--strand', dest='strand', default='+', help='strand to investigate matches, + or -')
    parser.add_argument('--embedded_txt', dest='embedded_txt', type=argparse.FileType('w'),
    default='embedded_results.txt', help='Results of potential embedded spanins')
    parser.add_argument('--overlap_txt', dest='overlap_txt', type=argparse.FileType('w'),
    default='overlap_results.txt', help='Results of potential overlapping spanins')
    parser.add_argument('--separate_txt', dest='separate_txt', type=argparse.FileType('w'),
    default='separated_results.txt', help='Results of potential separated spanins')
    parser.add_argument('--summary_txt', dest='summary_txt', type=argparse.FileType('w'),
    default='findSpanin_summary.txt', help='Results of potential spanin pairs')
    parser.add_argument('-v', action='version', version='0.3.0') # Is this manually updated?
    args = parser.parse_args()

    isp = getDescriptions(args.putative_isp_fasta_file)
    #print(len(isp))
    osp = getDescriptions(args.putative_osp_fasta_file)
    
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
    
    ###### The above steps probablt __SHOULD__ be wrapped into a little function. But, not necessary atm.

    # constructs list where we must multiply user input of AA by 3 to correspond to triplet codons
    embedded, overlap, separate = spaninProximity(data_isp,
                                                    data_osp, 
                                                    max_dist=args.max_isp_osp_distance*3) 
    s = 0
    for v in embedded.values():
        s += len(v)
    amt_embedded = s
    amt_unique_embedded = len(embedded.keys()) 
    s = 0
    for v in overlap.values():
        s += len(v)
    amt_overlap = s
    amt_unique_overlap = len(overlap.keys())
    s = 0
    for v in separate.values():
        s += len(v)
    amt_separate = s
    amt_unique_separate = len(separate.keys())

################################### OUTPUTS #################################################

    with args.summary_txt as f:
        f.write('++++++++++ Embedded Spanin Candidate Statistics +++++++++\n')
        f.writelines('Total Candidates = '+str(amt_embedded)+'\n')
        f.writelines('Unique ORF i-spanin = '+str(amt_unique_embedded))
        f.write('\n++++++++++ Overlap Spanin Candidate Statistics +++++++++\n')
        f.writelines('Total Candidates = '+str(amt_overlap)+'\n')
        f.writelines('Unique ORF i-spanin = '+str(amt_unique_overlap))
        f.write('\n++++++++++ Separate Spanin Candidate Statistics +++++++++\n')
        f.writelines('Total Candidates = '+str(amt_separate)+'\n')
        f.writelines('Unique ORF i-spanin = '+str(amt_unique_separate))
        f.write('\n++++++++++++++++++++++++ Totals +++++++++++++++++++++++++\n')
        f.writelines('Total Candidates = '+str(amt_embedded+amt_overlap+amt_separate))
        #f.writeline('Unique ORF i-spanin = '+str(amt_unique_embedded))

    with args.embedded_txt as f:
        f.write('================ Embedded Spanin Candidates =================\n')
        f.write('isp\tisp_start\tisp_end\tosp\tosp_start\tosp_end')
        if embedded != {}:
            for pisp, posp in embedded.items():
                f.write(pisp+'\n')
                for each_posp in posp:
                    f.write('\t{}\t{}\t{}\t{}\t{}\n'.format(each_posp[0],each_posp[1],each_posp[2],each_posp[3],each_posp[4]))
        else:
            f.write('nothing found')
    with args.overlap_txt as f:
        f.write('================ Overlap Spanin Candidates =================\n')
        f.write('isp\tisp_start\tisp_end\tosp\tosp_start\tosp_end')
        if overlap != {}:
            for pisp, posp in overlap.items():
                f.write(pisp+'\n')
                for each_posp in posp:
                    f.write('\t{}\t{}\t{}\t{}\t{}\n'.format(each_posp[0],each_posp[1],each_posp[2],each_posp[3],each_posp[4]))
        else:
            f.write('nothing found')

    with args.separate_txt as f:
        f.write('================ Separate Spanin Candidates =================\n')
        f.write('isp\tisp_start\tisp_end\tosp\tosp_start\tosp_end')
        if separate != {}:
            for pisp, posp in separate.items():
                f.write(pisp+'\n')
                for each_posp in posp:
                    f.write('\t{}\t{}\t{}\t{}\t{}\n'.format(each_posp[0],each_posp[1],each_posp[2],each_posp[3],each_posp[4]))
        else:
            f.write('nothing found')
    
