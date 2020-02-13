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

def extract_features(genbankFiles = None, fastaFiles = None, upOut = None, downOut = None, genesOnly = False, cdsOnly = True, forward = 1, behind = 1, outProt = True, tTable = 11, fTable = 11):

    genList = []
    fastaList = []
    
    for fileX in genbankFiles:        
        opener = SeqIO.parse(fileX, "genbank")
        for openRec in opener: # To turn the generator into objects (or else we end up with a list of generators)
            genList.append(openRec)

    for fileX in fastaFiles:
        opener = SeqIO.parse(fileX, "fasta")
        for openRec in opener: # Technically flattens multifastas too
            fastaList.append(openRec)


    for seqMatch in fastaList:
        longOut = seqMatch.description
        protID = seqMatch.id
        if(fTable != 0):
          fSeq = seqMatch.seq.translate(table = fTable, cds=False)
        else:
          fSeq = seqMatch.seq
        
        for gbk in genList:
            sourceOut = gbk.id
            num = -1
            for feat in gbk.features:
                num += 1
            
                if (genesOnly and feat.type != 'gene') or (cdsOnly and feat.type != 'CDS'):
                    continue

                temp = gbk.seq[feat.location.start : feat.location.end]
                if feat.location.strand == -1:
                    temp = temp.reverse_complement()

                if tTable != 0:
                  try:
                    gSeq = temp.translate(table = tTable, cds = True)
                  except CodonTable.TranslationError as cte:
                    #log.info("Translation issue at %s", cte)
                    gSeq = temp.translate(table=tTable, cds=False)
                else:
                  gSeq = temp

                if not('protein_id' in feat.qualifiers):
                    feat.qualifiers['protein_id'] = ["++++++++"] # Junk value for genesOnly flag

                if (gSeq == fSeq) or (protID == feat.qualifiers['protein_id'][0]): 
                    goBack = num - 1
                    goAhead = num + 1
                    numBack = behind
                    numAhead = forward
                    backList = []
                    aheadList = []

                    while (numBack != 0 and goBack >= 0):
                        if (genesOnly and gbk.features[goBack].type != 'gene') or (cdsOnly and gbk.features[goBack].type != 'CDS'):
                            goBack -= 1
                            continue
                        backList.append(gbk.features[goBack])
                        numBack -= 1
                        goBack -= 1

                    while (numAhead != 0 and goAhead < len(gbk.features)):
                        if (genesOnly and gbk.features[goAhead].type != 'gene') or (cdsOnly and gbk.features[goAhead].type != 'CDS'):
                            goAhead += 1
                            continue
                        aheadList.append(gbk.features[goAhead])
                        numAhead -= 1
                        goAhead += 1

                    
                    backList.reverse()
                
                    for item in backList:
                        if 'protein_id' in item.qualifiers:
                            upOut.write(">" + (item.qualifiers['protein_id'][0]) + " (5' of " + longOut + ' found within ' + sourceOut + ')\n')
                        else:
                            upOut.write(">" + (item.qualifiers['locus_tag'][0]) + " (5' of " + longOut + ' found within ' + sourceOut + ')\n')

                        if outProt == True:
                            if 'translation' in item.qualifiers:
                                upOut.write(str(item.qualifiers['translation'][0]) + '\n\n')
                            else:
                                modS = 0
                                modE = 0
                                if 'codon_start' in feature.qualifiers:
                                  if strand > 0:
                                    modS = int(feature.qualifiers['codon_start'][0]) - 1
                                  else:
                                    modE = int(feature.qualifiers['codon_start'][0]) - 1
                                seqHold = gbk.seq[item.location.start + modS: item.location.end + modE]
                                if item.location.strand == -1:
                                  seqHold = seqHold.reverse_complement()
                                if cdsOnly:
                                  if tTable != 0:
                                    upOut.write(str(seqHold.translate(table=tTable, cds=True)) + '\n\n')
                                  else:
                                    upOut.write(str(seqHold) + '\n\n')
                                else:
                                  if tTable != 0:
                                    upOut.write(str(seqHold.translate(table=tTable, cds=False)) + '\n\n')
                                  else:
                                    upOut.write(str(seqHold) + '\n\n')
                        else:
                            upOut.write(str(gbk.seq[item.location.start : item.location.end]) + '\n\n')

                    for item in aheadList:
                        if 'protein_id' in item.qualifiers:
                            downOut.write(">" + (item.qualifiers['protein_id'][0]) + " (3' of " + longOut + ' found within ' + sourceOut + ')\n')
                        else:
                            downOut.write(">" + (item.qualifiers['locus_tag'][0]) + " (3' of " + longOut + ' found within ' + sourceOut + ')\n')

                        if outProt == True:
                            if 'translation' in item.qualifiers:
                                downOut.write(str(item.qualifiers['translation'][0]) + '\n\n')
                            else:
                                modS = 0
                                modE = 0
                                if 'codon_start' in feature.qualifiers:
                                  if strand > 0:
                                    modS = int(feature.qualifiers['codon_start'][0]) - 1
                                  else:
                                    modE = int(feature.qualifiers['codon_start'][0]) - 1
                                seqHold = gbk.seq[item.location.start + modS: item.location.end + modE]
                                if item.location.strand == -1:
                                  seqHold = seqHold.reverse_complement()
                                if cdsOnly:
                                  if tTable != 0:
                                    downOut.write(str(seqHold.translate(table=tTable, cds=True)) + '\n\n')
                                  else:
                                    downOut.write(str(seqHold) + '\n\n')
                                else:
                                  if tTable != 0:
                                    downOut.write(str(seqHold.translate(table=tTable, cds=False)) + '\n\n')
                                  else:
                                    downOut.write(str(seqHold) + '\n\n')
                        else:
                            downOut.write(str(gbk.seq[item.location.start : item.location.end]) + '\n\n')
            #print(longOut)
    
    return

    
                    
 
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Export a subset of features from a Genbank file', epilog="")
    parser.add_argument('-genbankFiles', nargs='+', type=argparse.FileType("r"), help='Genbank file')
    parser.add_argument('-fastaFiles', nargs='+', type=argparse.FileType("r"), help='Fasta file to match against')
    parser.add_argument('-tTable', type=int, default=11, help='Translation table to use', choices=range(0, 23))
    parser.add_argument('-fTable', type=int, default=11, help='Translation table to use', choices=range(0, 23))
    parser.add_argument('-upOut', type=argparse.FileType("w"), help='Upstream Fasta output', default='test-data/upOut.fa')
    parser.add_argument('-downOut', type=argparse.FileType("w"), help='Downstream Fasta output', default='test-data/downOut.fa')
    parser.add_argument('--genesOnly', action='store_true', help='Search and return only Gene type features')
    parser.add_argument('--cdsOnly', action='store_true',  help='Search and return only CDS type features')
    parser.add_argument('--outProt', action='store_true',  help='Output the translated sequence')
    parser.add_argument('--forward', type=int, default=1, help='Number of features upstream from the hit to return')
    parser.add_argument('--behind', type=int, default=1, help='Number of features downstream from the hit to return')
    args = parser.parse_args()
    extract_features(**vars(args))
    
