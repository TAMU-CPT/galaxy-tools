#!/usr/bin/env python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_dna, generic_protein
import argparse
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

def extract_features(genbankFile = None, fastaFile = None, upOut = None, downOut = None, genesOnly = False, cdsOnly = True, numAhead = 2, numBack = 2):

    record = SeqIO.parse(genbankFile, "genbank")
    recordG = list(SeqIO.parse(fastaFile, "fasta"))
    matchF = recordG[0].seq

    for x in record:
        num = -1
        for y in x.features:
            num += 1
            
            if (genesOnly and y.type != 'gene') or (cdsOnly and y.type != 'CDS'):
                continue
            temp = x.seq[y.location.start : y.location.end]
            if y.location.strand == -1:
                temp = temp.reverse_complement()
            
            try:
                matchG = temp.translate(table = 11, cds = True)
            except CodonTable.TranslationError as cte:
                log.info("Translation issue at %s", cte)
                matchG = temp.translate(table=11, cds=False)
            if matchG == matchF:
                goBack = num - 1
                goAhead = num + 1
                backList = []
                aheadList = []
                while (numBack != 0 and goBack >= 0):
                    if (genesOnly and x.features[goBack].type != 'gene') or (cdsOnly and x.features[goBack].type != 'CDS'):
                        goBack -= 1
                        continue
                    backList.append(x.features[goBack])
                    numBack -= 1
                    goBack -= 1
                
                while (numAhead != 0 and goAhead < len(x.features)):
                    if (genesOnly and x.features[goAhead].type != 'gene') or (cdsOnly and x.features[goAhead].type != 'CDS'):
                        goAhead += 1
                        continue
                    aheadList.append(x.features[goAhead])
                    numAhead -= 1
                    goAhead += 1
                
                backList.reverse()
                
                for item in backList:
                    downOut.write(">" + (item.qualifiers['locus_tag'][0]) + '\n')
                    downOut.write(str(x.seq[item.location.start : item.location.end]) + '\n\n')
                
                for item in aheadList:
                    upOut.write(">" + (item.qualifiers['locus_tag'][0]) + '\n')
                    upOut.write(str(x.seq[item.location.start : item.location.end]) + '\n\n')
                    
 
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Export a subset of features from a Genbank file', epilog="")
    parser.add_argument('genbankFile', type=argparse.FileType("r"), help='Genbank file')
    parser.add_argument('fastaFile', type=argparse.FileType("r"), help='Fasta file to match against')
    parser.add_argument('upOut', type=argparse.FileType("w"), help='Upstream Fasta output', default='test-data/upOut.fa')
    parser.add_argument('downOut', type=argparse.FileType("w"), help='Downstream Fasta output', default='test-data/downOut.fa')
    parser.add_argument('--genesOnly', action='store_true', help='Search and return only Gene type features')
    parser.add_argument('--cdsOnly', action='store_true',  help='Search and return only CDS type features')
    parser.add_argument('--numAhead', type=int, default=2, help='Number of features upstream from the hit to return')
    parser.add_argument('--numBack', type=int, default=2, help='Number of features downstream from the hit to return')
    args = parser.parse_args()
    extract_features(**vars(args))
    
