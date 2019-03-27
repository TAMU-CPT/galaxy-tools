#!/usr/bin/env python

from BCBio import GFF
from Bio.SeqFeature import FeatureLocation
from Bio import SeqIO
from Bio.Seq import Seq
from gff3 import feature_lambda, feature_test_true
import csv
import argparse
import os #REMOVE BEFORE PUSH


def table_annotations(gff3In, tabularIn, fastaIn, out_gff3):

    if os.path.exists(os.getcwd()+"/test-data/changelog.txt"):
      os.remove(os.getcwd()+"/test-data/changelog.txt")
    if os.path.exists(os.getcwd()+"/test-data/snoke-edit.gff3"):
      os.remove(os.getcwd()+"/test-data/snoke-edit.gff3")

    # CSV parse tabular
    idDict = csv.DictReader(tabularIn, delimiter = '\t', fieldnames = ['OrgID', 'ID', 'Name', 'BoundS', 'BoundE', 'Strand', 'Note'])

    out_handle = open("test-data/snoke-edit.gff3", "a")
    # BioPython parse GFF
    
    #seq_dict = SeqIO.to_dict(SeqIO.parse(fastaIn, "fasta"))
    recG = list(GFF.parse(gff3In, SeqIO.to_dict(SeqIO.parse(fastaIn, "fasta"))))[0]

    # Get a changelog ready
    changelog = open("test-data/changelog.txt", "a")

    anyChange = False

    for row in idDict:
      Found = False;

      for i in recG.features:              
        if row["ID"] == i.id:
       
          if row["Strand"] == '+':
            row["Strand"] = +1
          else:
            row["Strand"] = -1

          Found = True
          Attr = False
   
          if i.qualifiers["Name"] != row["Name"]:
            i.qualifiers["Name"] = row["Name"]  
            Attr = True
         
          if i.location.start != int(row["BoundS"]) or i.location.end != int(row["BoundE"]) or i.strand != row["Strand"]:
            i.location = FeatureLocation(int(row["BoundS"]), int(row["BoundE"]), strand = row["Strand"] )
            Attr = True

          if "Note" in i.qualifiers and row["Note"] != "": # Set to not override with blank notes currently
            row["Note"] = (row["Note"]).split(',')
            if i.qualifiers["Note"] != row["Note"]:
              i.qualifiers["Note"] = row["Note"]
              Attr = True
          

          # On success, write out replaced attributes and success
          if Found and Attr:
            changelog.write(("Successfully overrode %s with \n" % i.id) + str(row) + '\n')
            anyChange = True

          # On fail, write out table line and why
          elif Found:
            # No changes detected
            changelog.write(("Did not override %d : No changes from table.\n Table line was " % i.id) + str(row) + '\n') 
      
      if Found == False:
        # No such ID
        changelog.write("Did not override %s : No such feature found.\n" % row["ID"])

    if anyChange:
      #out_handle = open("test-data/snoke-edit.gff3", "a")
      GFF.write([recG], out_gff3)
    else:
      changelog.write("No changes made to source GFF3 file.\n")

    changelog.close()
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Update GFF3 input from given tabular file')
    parser.add_argument('gff3In', type=argparse.FileType("r"), help='GFF3 source file')
    parser.add_argument('tabularIn', type=argparse.FileType("r"), help='Tabular file to fetch edits')
    parser.add_argument('fastaIn', type=argparse.FileType("r"), help='Tabular file to fetch edits')
    parser.add_argument('--out_gff3', type=argparse.FileType("w"), help='Output gff3', default='out.gff3')
    args = parser.parse_args()
    table_annotations(**vars(args))
