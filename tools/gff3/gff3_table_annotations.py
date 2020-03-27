#!/usr/bin/env python

from BCBio import GFF
from Bio.SeqFeature import FeatureLocation
from Bio import SeqIO
from Bio.Seq import Seq
import csv
import argparse


def table_annotations(gff3In, tabularIn, fastaIn, out_gff3, out_changelog):

    # CSV parse tabular
    header = csv.DictReader(tabularIn, delimiter="\t")
    for i in header.fieldnames:
        if i == "Boundary" and not ("bound_s" in header.fieldnames):
            header.fieldnames[header.fieldnames.index(i)] = "bound_s"

        elif i == "Boundary":
            header.fieldnames[header.fieldnames.index(i)] = "bound_e"
        # Else error
        elif i == "# Organism ID":
            header.fieldnames[header.fieldnames.index(i)] = "org_id"

        elif i == "User entered Notes":
            header.fieldnames[header.fieldnames.index(i)] = "Note"
    idDict = csv.DictReader(tabularIn, delimiter="\t", fieldnames=header.fieldnames)

    # BioPython parse GFF
    sourceG = list(GFF.parse(gff3In, SeqIO.to_dict(SeqIO.parse(fastaIn, "fasta"))))
    recG = []

    recTest = []

    sumFeatures = 0
    numOrgs = 0
    while numOrgs < len(sourceG):  # Should be directly editable
        topFeat = 0
        while topFeat < len(sourceG[numOrgs].features):
            subFeat1 = 0
            recG.append(sourceG[numOrgs].features[topFeat])
            sumFeatures += 1
            while subFeat1 < len(sourceG[numOrgs].features[topFeat].sub_features):
                subFeat2 = 0
                recG.append(sourceG[numOrgs].features[topFeat].sub_features[subFeat1])
                sumFeatures += 1
                while subFeat2 < len(
                    sourceG[numOrgs]
                    .features[topFeat]
                    .sub_features[subFeat1]
                    .sub_features
                ):
                    recG.append(
                        sourceG[numOrgs]
                        .features[topFeat]
                        .sub_features[subFeat1]
                        .sub_features[subFeat2]
                    )
                    sumFeatures += 1
                    subFeat2 += 1
                subFeat1 += 1
            topFeat += 1
        numOrgs += 1

    # Get a changelog ready
    out_changelog.write("ID\tChanges\tStatus\n")

    anyChange = False

    for row in idDict:
        if row["ID"] == "ID":
            continue  # Skip header
        Found = False

        for i in recG:

            # if "Parent" in i.qualifiers:

            if row["ID"] == i.id:

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
                    if qual == "ID" or qual == "":
                        continue

                    if qual == "Strand":
                        if row["Strand"] == "+":
                            row["Strand"] = +1
                        else:
                            row["Strand"] = -1

                    if qual == "Name" and row["Name"] != i.qualifiers["Name"][0]:
                        i.qualifiers["Name"][0] = row["Name"]
                        nameC = True

                    elif qual == "Alias" and row["Alias"] != i.qualifiers["Alias"][0]:
                        i.qualifiers["Alias"][0] = row["Alias"]
                        aliasC = True

                    elif (
                        qual == "Parent" and row["Parent"] != i.qualifiers["Parent"][0]
                    ):
                        i.qualifiers["Parent"][0] = row["Parent"]
                        parentC = true

                    # elif qual == "Note":

                    elif qual == "Strand" and i.strand != row["Strand"]:
                        strandC = True
                        i.location = FeatureLocation(
                            i.location.start, i.location.end, row["Strand"]
                        )

                    elif qual == "bound_s" and i.location.start != int(row["bound_s"]):
                        startC = True
                        i.location = FeatureLocation(
                            int(row["bound_s"]), i.location.end, i.location.strand
                        )

                    elif qual == "bound_e" and i.location.end != int(row["bound_e"]):
                        endC = True
                        i.location = FeatureLocation(
                            i.location.start, int(row["bound_e"]), i.location.strand
                        )

                    elif qual == "Note":
                        temp = [str(row["Note"])]
                        if ("Note" in i.qualifiers) and i.qualifiers["Note"] != temp:
                            if temp:
                                i.qualifiers["Note"] = temp
                            else:
                                i.qualifiers.pop("Note", None)
                            noteC = True
                        elif temp != [""] and not ("Note" in i.qualifiers):
                            i.qualifiers["Note"] = temp
                            noteC = True

                    elif not (
                        qual
                        in [
                            "Target",
                            "Gap",
                            "Dbxref",
                            "Ontology_term",
                            "Is_circular",
                            "Derives_from",
                            "bound_s",
                            "bound_e",
                            "org_id",
                            "Strand",
                            "Name",
                            "Note",
                        ]
                    ):
                        temp = qual.lower().replace(" ", "_")
                        if temp in i.qualifiers:  # Edit
                            if type(row[qual]) == type(None):
                                i.qualifiers.pop(temp, None)
                                qualC = True
                            elif i.qualifiers[temp] != [str(row[qual])]:
                                i.qualifiers[temp] = [str(row[qual])]
                                qualC = True
                        elif type(row[qual]) != type(None):  # Create
                            i.qualifiers[temp] = [str(row[qual])]
                            qualC = True

                    # print(i)
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

                if (
                    changeList != ""
                ):  # On success, write out replaced attributes and success
                    out_changelog.write("%s\t%s\tSuccess\n" % (i.id, changeList))
                    anyChange = True
                else:  # On fail, write out table line and why
                    # No changes detected
                    out_changelog.write("%s\tNone\tNo Change\n" % i.id)

                break

        if Found == False:
            # No such ID
            out_changelog.write("%s\tNone\tID not Found\n" % row["ID"])

    if anyChange:
        sourceG[0].annotations = {}
        sourceG[0].features = [x for x in sourceG[0].features if x.type != "remark"]
        GFF.write(sourceG, out_gff3)
    else:
        out_changelog.write("GFF3\tNone\tGFF3 already equals Table\n")
        out_gff3 = gff3In

    out_changelog.close()
    out_gff3.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Update GFF3 input from given tabular file"
    )
    parser.add_argument("gff3In", type=argparse.FileType("r"), help="GFF3 source file")
    parser.add_argument(
        "tabularIn", type=argparse.FileType("r"), help="Tabular file to fetch edits"
    )
    parser.add_argument(
        "fastaIn",
        type=argparse.FileType("r"),
        help="Fasta file to correctly construct the GFF3",
    )
    parser.add_argument(
        "--out_gff3",
        type=argparse.FileType("w"),
        help="Output gff3",
        default="out.gff3",
    )
    parser.add_argument(
        "--out_changelog",
        type=argparse.FileType("w"),
        help="Output changelog",
        default="out.tsv",
    )
    args = parser.parse_args()
    table_annotations(**vars(args))  # 283
