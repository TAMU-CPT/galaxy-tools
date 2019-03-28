#!/usr/bin/env python

from BCBio import GFF
from Bio.SeqFeature import FeatureLocation
from Bio import SeqIO
from Bio.Seq import Seq
from gff3 import feature_lambda, feature_test_true
import csv
import argparse

def table_annotations(gff3In, tabularIn, fastaIn, out_gff3, out_changelog):

    # CSV parse tabular
    idDict = csv.DictReader(tabularIn, delimiter = '\t', fieldnames = ['OrgID', 'ID', 'Name', 'BoundS', 'BoundE', 'Strand', 'Note'])

    # BioPython parse GFF
    recG = list(GFF.parse(gff3In, SeqIO.to_dict(SeqIO.parse(fastaIn, "fasta"))))[0]

    # Get a changelog ready
    out_changelog.write("ID\tChanges\tStatus\n")

    anyChange = False

    for row in idDict:
      if row["ID"] == "ID":
        continue # Skip header
      Found = False;

      for i in recG.features:              
        if row["ID"] == i.id:
          
          strandC = False
          startC = False
          endC = False
          nameC = False
          noteC = False
          
          if row["Strand"] == '+':
            row["Strand"] = +1
          else:
            row["Strand"] = -1

          Found = True
          Attr = False

          if (row["Name"] != i.qualifiers["Name"][0]):
            i.qualifiers["Name"][0] = row["Name"]  
            nameC = True
         
          if i.location.start != int(row["BoundS"]):
            startC = True
          
          if i.location.end != int(row["BoundE"]):
            endC = True
 
          if i.strand != row["Strand"]:
            strandC = True

          if startC or endC or strandC: # Location object needs to be rebuilt, can't individually set start/end
            i.location = FeatureLocation(int(row["BoundS"]), int(row["BoundE"]), strand = row["Strand"] )

          if ("Note" in i.qualifiers) and row["Note"] != "": # Set to not override with blank notes currently
            row["Note"] = (row["Note"]).split(',') # Turn note into a list
            if i.qualifiers["Note"] != row["Note"]:
              i.qualifiers["Note"] = row["Note"]
              noteC = True
          elif row["Note"] != "":
            row["Note"] = (row["Note"]).split(',')
            i.qualifiers["Note"] = row["Note"]
            noteC = True
          
          changeList = ""
          if nameC:
            changeList += "Name"
          if startC:
            if changeList != "":
              changeList += ", "
            changeList += "Start"            
          if endC:
            if changeList != "":
              changeList += ", "
            changeList += "End" 
          if strandC:
            if changeList != "":
              changeList += ", "
            changeList += "Strand" 
          if noteC:
            if changeList != "":
              changeList += ", "
            changeList += "Notes" 
          
          if changeList != "": # On success, write out replaced attributes and success
            out_changelog.write("%s\t%s\tSuccess\n" % (i.id, changeList))
            anyChange = True          
          else: # On fail, write out table line and why
            # No changes detected
            out_changelog.write("%s\tNone\tNo Change\n" % i.id) 
      
      if Found == False:
        # No such ID
        out_changelog.write("%s\tNone\tID not Found\n" % row["ID"])

    if anyChange:
      #out_handle = open("test-data/snoke-edit.gff3", "a")
      recG.features = [x for x in recG.features if x.type != 'remark']
      GFF.write([recG], out_gff3)
    else:
      out_changelog.write("GFF3\tNone\tGFF3 already equals Table\n")

    out_changelog.close()
    out_gff3.close()
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Update GFF3 input from given tabular file')
    parser.add_argument('gff3In', type=argparse.FileType("r"), help='GFF3 source file')
    parser.add_argument('tabularIn', type=argparse.FileType("r"), help='Tabular file to fetch edits')
    parser.add_argument('fastaIn', type=argparse.FileType("r"), help='Fasta file to correctly construct the GFF3')
    parser.add_argument('--out_gff3', type=argparse.FileType("w"), help='Output gff3', default='out.gff3')
    parser.add_argument('--out_changelog', type=argparse.FileType("w"), help='Output changelog', default='out.tsv')
    args = parser.parse_args()
    table_annotations(**vars(args))
