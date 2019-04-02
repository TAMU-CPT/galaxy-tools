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
    header = csv.DictReader(tabularIn, delimiter = '\t')
    for i in header.fieldnames:
      if(i == "Boundary" and not("bound_s" in header.fieldnames)):
        header.fieldnames[header.fieldnames.index(i)] = "bound_s"

      elif(i == "Boundary"):
        header.fieldnames[header.fieldnames.index(i)] = "bound_e"
      #Else error
      elif(i == "# Organism ID"):
        header.fieldnames[header.fieldnames.index(i)] = "org_id"

      elif(i == "User entered Notes"):
        header.fieldnames[header.fieldnames.index(i)] = "Note"
    idDict = csv.DictReader(tabularIn, delimiter = '\t', fieldnames = header.fieldnames)

    # BioPython parse GFF
    sourceG = list(GFF.parse(gff3In, SeqIO.to_dict(SeqIO.parse(fastaIn, "fasta"))))
    recG = []   
    
    parentDict = {}

    topCount = 0
    
    for i in sourceG:
      topCount += 1
      if i.features:
        recG.extend(i.features) # Adds Genes
        for j in i.features: # For each gene
          sub1 = j.sub_features #sub1 = mRNA
          if sub1:
            parentDict[j.id] = j
          recG.extend(sub1)
          for k in j.sub_features:
            if k.sub_features: 
              parentDict[k.id] = k
              recG.extend(k.sub_features)            
 
    # Get a changelog ready
    out_changelog.write("ID\tChanges\tStatus\n")

    anyChange = False

    for row in idDict:
      if row["ID"] == "ID":
        continue # Skip header
      Found = False
        

      for i in recG:

        #if "Parent" in i.qualifiers:
        

        if row["ID"] == i.id:

          numLevels = 0
          treeTraverse = ["", "", ""]     
          tree = i     
          topGene = ""
          while "Parent" in tree.qualifiers:
            treeTraverse[numLevels] = tree 
            numLevels += 1
            tree = parentDict[tree.qualifiers["Parent"]]
            
          treeTraverse[numLevels] = tree
 
          change = sourceG[0]

          if numLevels == 2:
            for find in change.features:
              if treeTraverse[2].id == find.id:
                change = find
                break
            for find in change.sub_features:
              if treeTraverse[1].id == find.id:
                change = find
                break
            for find in change.sub_features:
              if treeTraverse[0].id == find.id:
                change = find
                break

          elif numLevels == 1:
            for find in change.features:
              if treeTraverse[1].id == find.id:
                change = find
                break
            for find in change.sub_features:
              if treeTraverse[0].id == find.id:
                change = find 
                break     

          if numLevels == 0:
            for find in change.features:
              if treeTraverse[0].id == find.id:
                change = find
                break

          strandC = False
          startC = False
          endC = False
          nameC = False
          noteC = False
          qualC = False
          aliasC = False
          parentC = False

          Found = True

          

          for qual in row:
            if qual == "ID":
              continue

            if qual == "Strand":
              if row["Strand"] == '+':
                row["Strand"] = +1
              else:
                row["Strand"] = -1
            
            if qual == "Name" and row["Name"] != change.qualifiers["Name"][0]:
              change.qualifiers["Name"][0] = row["Name"]
              nameC = True

            elif qual == "Alias" and row["Alias"] != change.qualifiers["Alias"]:
              change.qualifiers["Alias"][0] = row["Alias"]
              aliasC = True

            elif qual == "Parent" and row["Parent"] != change.qualifiers["Parent"][0]:
              change.qualifiers["Parent"][0] = row["Parent"]
              parentC = true

           # elif qual == "Note":

            elif qual == "Strand" and change.strand != row["Strand"]:
              strandC = True
              change.location = FeatureLocation(change.location.start, change.location.end, row["Strand"])
          

            elif qual == "bound_s" and change.location.start != int(row["bound_s"]):
              startC = True
              change.location = FeatureLocation(int(row["bound_s"]), change.location.end, change.location.strand)

            elif qual == "bound_e" and change.location.end != int(row["bound_e"]):
              endC = True
              change.location = FeatureLocation(change.location.start, int(row["bound_e"]), change.location.strand)
            
            elif qual == "Note":
              temp = [str(row["Note"])]
              if ("Note" in change.qualifiers) and change.qualifiers["Note"] != temp:
                if temp:
                  change.qualifiers["Note"] = temp
                else:
                  change.qualifiers.pop("Note", None)
                noteC = True
              elif temp != [''] and not ("Note" in change.qualifiers):
                change.qualifiers["Note"] = temp
                noteC = True

            elif not (qual in ['Target', 'Gap', 'Dbxref', 'Ontology_term', 'Is_circular', 'Derives_from', 'bound_s', 'bound_e', 'org_id', 'Strand', 'Name', 'Note']):
              temp = qual.lower().replace(' ', '_')
              if temp in change.qualifiers: # Edit
                if type(row[qual]) == type(None):
                  change.qualifiers.pop(temp, None)
                  qualC = True
                elif change.qualifiers[temp] != [str(row[qual])]:
                  change.qualifiers[temp] = [str(row[qual])]
                  qualC = True
              elif type(row[qual]) != type(None): # Create
                change.qualifiers[temp] = [str(row[qual])]
                qualC = True
            
            #print(i)
          # if "OrgID" in row and (row["OrgID"] != i.qualifiers["Name"][0]):
          # OrgID Seemingly not used aside from GFF Header

          # Location object needs to be rebuilt, can't individually set start/end

          changeList = ""
          if nameC:
            changeList += "Name"
          if aliasC:
            if changeList != "":
              changeList += ", "
            changeList += "Alias"
          if parentC:
            if changeList != "":
              changeList += ", "
            changeList += "Parent"
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
          if qualC:
            if changeList != "":
              changeList += ", "
            changeList += "Other Qualifiers"

          if changeList != "": # On success, write out replaced attributes and success
            out_changelog.write("%s\t%s\tSuccess\n" % (change.id, changeList))
            anyChange = True
          else: # On fail, write out table line and why
            # No changes detected
            out_changelog.write("%s\tNone\tNo Change\n" % change.id)

          break

      if Found == False:
        # No such ID
        out_changelog.write("%s\tNone\tID not Found\n" % row["ID"])

    if anyChange:
      sourceG[0].annotations = {}
      sourceG[0].features = [x for x in sourceG[0].features if x.type != 'remark']
      GFF.write(sourceG, out_gff3)
    else:
      out_changelog.write("GFF3\tNone\tGFF3 already equals Table\n")
      out_gff3 = gff3In

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
