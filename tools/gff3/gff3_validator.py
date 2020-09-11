#!/usr/bin/env python

from cpt_gffParser import gffParse, gffWrite
from Bio.SeqFeature import FeatureLocation
from Bio import SeqIO, SeqFeature, SeqRecord
from Bio.Seq import Seq
from gff3 import feature_lambda, feature_test_true
import csv
import argparse
from cpt_gffParser import gffParse, gffWrite

def table_annotations(gff3In, out_errorlog):

    numError = 0
    numWarning = 0
    containsHash = {}
    locations = []
    # disallowedChars = ['\t', '\n', '\r', '%', ';', '=', '&',',']
    disallowedChars = ['"']
    #"""
    outList = (gffParse(gff3In))
    gffWrite(None)
    #GFFWrite(outList)
    #for x in outList:
    #  GFFWrite(x)
      #if "sequence-region" in x.annotations.keys():
      #  print(x.annotations["sequence-region"][2])
      #for y in x.features:
      #  if y.sub_features:
      #    for z in y.sub_features:
      #      if z.sub_features:
      #        for t in z.sub_features:
      #          print(t)
    exit()
    #"""
    # For each record in gff3 parse
    for record in list(gffParse(gff3In)):
        print((record.annotations))
        # For each feature in record
        # featLvl1 = Genes, terminator, tRNA
        # featLvl2 = mRNA
        # featLvl3 = Exons, CDS, introns, and Shines
        #for x in record.features:
        #  print(x)
        exit()
        errorMessage = ""

        for featLvl1 in record.features:
            for notes in featLvl1.qualifiers:
                for i in disallowedChars:
                    for j in featLvl1.qualifiers[notes]:
                        if i in j:
                            if i == "\n":
                                problem = "New Line (Return)"
                            elif i == "\t":
                                problem = "Tab space"
                            elif i == "\r":
                                problem = "Carriage return"
                            else:
                                problem = i
                            # errorMessage = errorMessage + ("Warning: Unencoded character in qualifiers of %s %s\n  Problem qualifier: %s: %s --- Unencoded Character: %s\n\n" % (featLvl1.type, featLvl1.id, str(notes), j, problem))
                            errorMessage = errorMessage + (
                                "Warning: Character '\"' in qualifiers of %s %s, will not be able to convert to Genbank\n Problem qualifier: %s: %s\n\n"
                                % (featLvl1.type, featLvl1.id, str(notes), j)
                            )
                            numWarning += 1

            dupCheck = containsHash.get(
                (
                    featLvl1.type
                    + str(featLvl1.location.start)
                    + str(featLvl1.location.end)
                ),
                -1,
            )
            if dupCheck != -1:
                errorMessage = errorMessage + (
                    "Error: Duplicate %ss in GFF3 file\n  %s ID %s occupies the same space as %s ID: %s\n  Note: if one of these IDs doesn't appear to exist, ensure there are not two of the other\n\n"
                    % (
                        featLvl1.type,
                        dupCheck.type,
                        dupCheck.id,
                        featLvl1.type,
                        featLvl1.id,
                    )
                )
                numError += 1
            else:
                containsHash[
                    (
                        featLvl1.type
                        + str(featLvl1.location.start)
                        + str(featLvl1.location.end)
                    )
                ] = featLvl1

            foundMRNA = 0
            locations.append(featLvl1)

            # For each gene
            if featLvl1.type == "gene":

                for featLvl2 in featLvl1.sub_features:

                    for notes in featLvl2.qualifiers:
                        for i in disallowedChars:
                            for j in featLvl2.qualifiers[notes]:
                                if i in j:
                                    if i == "\n":
                                        problem = "New Line (Return)"
                                    elif i == "\t":
                                        problem = "Tab space"
                                    elif i == "\r":
                                        problem = "Carriage return"
                                    else:
                                        problem = i
                                    # errorMessage = errorMessage + ("Warning: Unencoded character in qualifiers of %s %s\n  Problem qualifier: %s: %s --- Unencoded Character: %s\n\n" % (featLvl2.type, featLvl2.id, str(notes), j, problem))
                                    errorMessage = errorMessage + (
                                        "Warning: Character '\"' in qualifiers of %s %s, will not be able to convert to Genbank\n Problem qualifier: %s: %s\n\n"
                                        % (featLvl1.type, featLvl1.id, str(notes), j)
                                    )
                                    numWarning += 1

                    dupCheck = containsHash.get(
                        (
                            featLvl2.type
                            + str(featLvl2.location.start)
                            + str(featLvl2.location.end)
                        ),
                        -1,
                    )
                    if dupCheck != -1:
                        errorMessage = errorMessage + (
                            "Error: Duplicate %ss in GFF3 file\n  %s ID %s occupies the same space as %s ID: %s\n  Note: if one of these IDs doesn't appear to exist, ensure there are not two of the other\n\n"
                            % (
                                featLvl2.type,
                                dupCheck.type,
                                dupCheck.id,
                                featLvl2.type,
                                featLvl2.id,
                            )
                        )
                        numError += 1
                    else:
                        containsHash[
                            (
                                featLvl2.type
                                + str(featLvl2.location.start)
                                + str(featLvl2.location.end)
                            )
                        ] = featLvl2

        numError = 0
        numWarning = 0
        containsHash = {}
        locations = []
        disallowedChars = ["\t", "\n", "\r", "%", ";", "=", "&", ","]

        # For each record in gff3 parse
        for record in list(gffParse(gff3In)):
            # For each feature in record
            # featLvl1 = Genes, terminator, tRNA
            # featLvl2 = mRNA
            # featLvl3 = Exons, CDS, introns, and Shines

            errorMessage = ""

            for featLvl1 in record.features:

                for notes in featLvl1.qualifiers:
                    for i in disallowedChars:
                        for j in featLvl3.qualifiers[notes]:
                            if i in j:
                                if i == "\n":
                                    problem = "New Line (Return)"
                                elif i == "\t":
                                    problem = "Tab space"
                                elif i == "\r":
                                    problem = "Carriage return"
                                else:
                                    problem = i
                                # errorMessage = errorMessage + ("Warning: Unencoded character in qualifiers of %s %s\n  Problem qualifier: %s: %s --- Unencoded Character: %s\n\n" % (featLvl3.type, featLvl3.id, str(notes), j, problem))
                                errorMessage = errorMessage + (
                                    "Warning: Character '\"' in qualifiers of %s %s, will not be able to convert to Genbank\n Problem qualifier: %s: %s\n\n"
                                    % (featLvl1.type, featLvl1.id, str(notes), j)
                                )
                                numWarning += 1

                    # cdsList = []
                    # exonList = []

                dupCheck = containsHash.get(
                    (
                        featLvl3.type
                        + str(featLvl3.location.start)
                        + str(featLvl3.location.end)
                    ),
                    -1,
                )
                if dupCheck != -1:
                    errorMessage = errorMessage + (
                        "Error: Duplicate %ss in GFF3 file\n  %s ID %s occupies the same space as %s ID: %s\n  Note: if one of these IDs doesn't appear to exist, ensure there are not two of the other\n\n"
                        % (
                            featLvl3.type,
                            dupCheck.type,
                            dupCheck.id,
                            featLvl3.type,
                            featLvl3.id,
                        )
                    )
                    numError += 1
                else:
                    containsHash[
                        (
                            featLvl1.type
                            + str(featLvl1.location.start)
                            + str(featLvl1.location.end)
                        )
                    ] = featLvl1

                foundMRNA = 0
                locations.append(featLvl1)

                # For each gene
                if featLvl1.type == "gene":

                    for featLvl2 in featLvl1.sub_features:

                        for notes in featLvl2.qualifiers:
                            for i in disallowedChars:
                                for j in featLvl2.qualifiers[notes]:
                                    if i in j:
                                        if i == "\n":
                                            problem = "New Line (Return)"
                                        elif i == "\t":
                                            problem = "Tab space"
                                        elif i == "\r":
                                            problem = "Carriage return"
                                        else:
                                            problem = i
                                        errorMessage = errorMessage + (
                                            "Warning: Unencoded character in qualifiers of %s %s\n  Problem qualifier: %s: %s --- Unencoded Character: %s\n\n"
                                            % (
                                                featLvl2.type,
                                                featLvl2.id,
                                                str(notes),
                                                j,
                                                problem,
                                            )
                                        )
                                        numWarning += 1

                        dupCheck = containsHash.get(
                            (
                                featLvl2.type
                                + str(featLvl2.location.start)
                                + str(featLvl2.location.end)
                            ),
                            -1,
                        )
                        if dupCheck != -1:
                            errorMessage = errorMessage + (
                                "Error: Duplicate %ss in GFF3 file\n  %s ID %s occupies the same space as %s ID: %s\n  Note: if one of these IDs doesn't appear to exist, ensure there are not two of the other\n\n"
                                % (
                                    featLvl2.type,
                                    dupCheck.type,
                                    dupCheck.id,
                                    featLvl2.type,
                                    featLvl2.id,
                                )
                            )
                            numError += 1
                        else:
                            containsHash[
                                (
                                    featLvl2.type
                                    + str(featLvl2.location.start)
                                    + str(featLvl2.location.end)
                                )
                            ] = featLvl2

                        if featLvl2.type == "mRNA":
                            foundMRNA += 1

                            if featLvl2.location.start != featLvl1.location.start:
                                errorMessage = (
                                    errorMessage
                                    + "Error: mRNA start boundary does not match parent gene's boundary\n  Gene %s Start Boundary: %d\n  mRNA %s Start Boundary: %d\n\n"
                                    % (
                                        featLvl1.id,
                                        featLvl1.location.start,
                                        featLvl2.id,
                                        featLvl2.location.start,
                                    )
                                )
                                numError += 1

                            if featLvl2.location.end != featLvl1.location.end:
                                errorMessage = (
                                    errorMessage
                                    + "Error: mRNA end boundary does not match parent gene's boundary\n  Gene %s End Boundary: %d\n  mRNA %s End Boundary: %d\n\n"
                                    % (
                                        featLvl1.id,
                                        featLvl1.location.end,
                                        featLvl2.id,
                                        featLvl2.location.end,
                                    )
                                )
                                numError += 1

                            exons = 0
                            CDS = 0
                            shines = 0
                            cdsList = []
                            exonList = []

                            for featLvl3 in featLvl2.sub_features:

                                for notes in featLvl3.qualifiers:
                                    for i in disallowedChars:
                                        for j in featLvl3.qualifiers[notes]:
                                            if i in j:
                                                if i == "\n":
                                                    problem = "New Line (Return)"
                                                elif i == "\t":
                                                    problem = "Tab space"
                                                elif i == "\r":
                                                    problem = "Carriage return"
                                                else:
                                                    problem = i
                                                errorMessage = errorMessage + (
                                                    "Warning: Unencoded character in qualifiers of %s %s\n  Problem qualifier: %s: %s --- Unencoded Character: %s\n\n"
                                                    % (
                                                        featLvl3.type,
                                                        featLvl3.id,
                                                        str(notes),
                                                        j,
                                                        problem,
                                                    )
                                                )
                                                numWarning += 1

                                # cdsList = []
                                # exonList = []

                                dupCheck = containsHash.get(
                                    (
                                        featLvl3.type
                                        + str(featLvl3.location.start)
                                        + str(featLvl3.location.end)
                                    ),
                                    -1,
                                )
                                if dupCheck != -1:
                                    errorMessage = errorMessage + (
                                        "Error: Duplicate %ss in GFF3 file\n  %s ID %s occupies the same space as %s ID: %s\n  Note: if one of these IDs doesn't appear to exist, ensure there are not two of the other\n\n"
                                        % (
                                            featLvl3.type,
                                            dupCheck.type,
                                            dupCheck.id,
                                            featLvl3.type,
                                            featLvl3.id,
                                        )
                                    )
                                    numError += 1
                                else:
                                    containsHash[
                                        (
                                            featLvl3.type
                                            + str(featLvl3.location.start)
                                            + str(featLvl3.location.end)
                                        )
                                    ] = featLvl3

                                if featLvl3.type == "exon":
                                    exons += 1
                                    exonList.append(featLvl3)
                                    if (
                                        featLvl3.location.start
                                        < featLvl2.location.start
                                        or featLvl3.location.start
                                        > featLvl2.location.end
                                    ):
                                        errorMessage += (
                                            "Error: Exon start falls outside the boundary of parent mRNA\n  %s mRNA Start: %d  -- End: %d\n  %s Exon Start: %d\n\n"
                                            % (
                                                featLvl2.id,
                                                featLvl2.location.start,
                                                featLvl2.location.end,
                                                featLvl3.id,
                                                featLvl3.location.start,
                                            )
                                        )
                                        numError += 1
                                    elif (
                                        featLvl3.location.end < featLvl2.location.start
                                        or featLvl3.location.end > featLvl2.location.end
                                    ):
                                        errorMessage += (
                                            "Error: Exon end falls outside the boundary of parent mRNA\n  %s mRNA Start: %d  -- End: %d\n  %s Exon End: %d\n\n"
                                            % (
                                                featLvl2.id,
                                                featLvl2.location.start,
                                                featLvl2.location.end,
                                                featLvl3.id,
                                                featLvl3.location.end,
                                            )
                                        )
                                        numError += 1
                                    elif (
                                        featLvl3.location.start
                                        != featLvl2.location.start
                                        and featLvl3.location.end
                                        != featLvl2.location.end
                                    ):
                                        errorMessage += (
                                            "Error: Exon does not match expected start/end pair based on parent mRNA boundary\n  %s mRNA Start: %d  -- End: %d\n  %s Exon Start: %d  -- End: %d\n\n"
                                            % (
                                                featLvl2.id,
                                                featLvl2.location.start,
                                                featLvl2.location.end,
                                                featLvl3.id,
                                                featLvl3.location.start,
                                                featLvl3.location.end,
                                            )
                                        )
                                        numError += 1

                                elif featLvl3.type == "CDS":
                                    CDS += 1
                                    cdsList.append(featLvl3)

                                    if (
                                        featLvl3.location.start
                                        < featLvl2.location.start
                                        or featLvl3.location.start
                                        > featLvl2.location.end
                                    ):
                                        errorMessage += (
                                            "Error: CDS start falls outside the boundary of parent mRNA\n  %s mRNA Start: %d  -- End: %d\n  %s CDS Start: %d\n\n"
                                            % (
                                                featLvl2.id,
                                                featLvl2.location.start,
                                                featLvl2.location.end,
                                                featLvl3.id,
                                                featLvl3.location.start,
                                            )
                                        )
                                        numError += 1
                                    elif (
                                        featLvl3.location.end < featLvl2.location.start
                                        or featLvl3.location.end > featLvl2.location.end
                                    ):
                                        errorMessage += (
                                            "Error: CDS end falls outside the boundary of parent mRNA\n  %s mRNA Start: %d  -- End: %d\n  %s CDS End: %d\n\n"
                                            % (
                                                featLvl2.id,
                                                featLvl2.location.start,
                                                featLvl2.location.end,
                                                featLvl3.id,
                                                featLvl3.location.end,
                                            )
                                        )
                                        numError += 1

                            #              elif(featLvl3.type == "Shine_Dalgarno_sequence"):
                            #                shines += 1 # SHINE GET!
                            #                if(featLvl3.location.start < featLvl2.location.start or featLvl3.location.start > featLvl2.location.end):
                            #                  errorMessage += ("Error: Shine Dalgarno start falls outside the boundary of parent mRNA\n  %s mRNA Start: %d  -- End: %d\n  %s Shine Dalgarno Start: %d\n\n" % (featLvl2.id, featLvl2.location.start, featLvl2.location.end, featLvl3.id, featLvl3.location.start))
                            #                elif(featLvl3.location.end < featLvl2.location.start or featLvl3.location.end > featLvl2.location.end):
                            #                  errorMessage += ("Error: Shine Dalgarno end falls outside the boundary of parent mRNA\n  %s mRNA Start: %d  -- End: %d\n  %s Shine Dalgarno End: %d\n\n" % (featLvl2.id, featLvl2.location.start, featLvl2.location.end, featLvl3.id, featLvl3.location.end))
                            #                elif(featLvl3.location.start != featLvl2.location.start and featLvl3.location.end != featLvl2.location.end):
                            #                  errorMessage += ("Error: Shine Dalgarno does not match expected start/end pair based on parent mRNA boundary\n  %s mRNA Start: %d  -- End: %d\n  %s Shine Dalgarno Start: %d  -- End: %d\n\n" % (featLvl2.id, featLvl2.location.start, featLvl2.location.end, featLvl3.id, featLvl3.location.start, featLvl3.location.end))

                            if (
                                cdsList
                            ):  # Interpreter requires this, can't just check against CDS sum
                                if CDS == 1:
                                    cdsCheck = cdsList.pop()
                                    if (
                                        cdsCheck.location.start
                                        != featLvl2.location.start
                                        and cdsCheck.location.end
                                        != featLvl2.location.end
                                    ):
                                        errorMessage += (
                                            "Error: CDS does not match expected start/end pair based on parent mRNA boundary\n  %s mRNA Start: %d  -- End: %d\n  %s CDS Start: %d  -- End: %d\n\n"
                                            % (
                                                featLvl2.id,
                                                featLvl2.location.start,
                                                featLvl2.location.end,
                                                cdsList[0].id,
                                                cdsList[0].location.start,
                                                cdsList[0].location.end,
                                            )
                                        )
                                        numError += 1
                                elif CDS > 1:
                                    match = False
                                    for i in CDS:
                                        cdsCheck = cdsList.pop()
                                        if (
                                            cdsCheck.location.start
                                            != featLvl2.location.start
                                            or cdsCheck.location.end
                                            != featLvl2.location.end
                                        ):
                                            match = True
                                    if match == False:
                                        errorMessage = (
                                            "Error: Multiple CDS in mRNA %s, but none match expected start/end boundary\n\n"
                                            % (featLvl2.id)
                                        ) + errorMessage
                                        numError += 1

                            else:
                                errorMessage = (
                                    "Warning: %s mRNA has no CDS features\n\n"
                                    % (featLvl2.id)
                                ) + errorMessage
                                numWarning += 1

                            if exons > 2:
                                errorMessage = (
                                    "Error: %s mRNA has %d exons, expected no more than 2\n\n"
                                    % (featLvl2.id, exons)
                                ) + errorMessage
                                numError += 1
                            elif exons > 1 and len(exonList) > 1:
                                if not (
                                    exonList[0].location.start
                                    > exonList[1].location.end
                                    or exonList[1].location.start
                                    > exonList[0].location.end
                                ):
                                    errorMessage = (
                                        "Error: mRNA has 2 exons, and they don't match the expected start/end boundaries\n  mRNA %s Start: %d -- End: %d\n  Exon %s Start: %d -- End: %d\n  Exon %s Start: %d -- End: %d\n  Expected one exon's start to match the mRNA's start, the other to match its end, and neither to overlap.\n\n"
                                        % (
                                            featLvl2.id,
                                            featLvl2.location.start,
                                            featLvl2.location.end,
                                            exonList[0].id,
                                            exonList[0].location.start,
                                            exonList[0].location.end,
                                            exonList[1].id,
                                            exonList[1].location.start,
                                            exonList[1].location.end,
                                        )
                                    ) + errorMessage
                                    numError += 1
                            elif exons == 0:
                                errorMessage = (
                                    "Warning: %s mRNA has no exons\n\n" % (featLvl2.id)
                                ) + errorMessage
                                numWarning += 1
                            #            if(shines > 0):
                            #              errorMessage = ("Error: %s mRNA has %d Shine Dalgarno sequences, expected no more than 1\n\n" % (featLvl2.id, shines))  + errorMessage

                            if exons + CDS + shines == 0:
                                errorMessage = (
                                    "Warning: %s mRNA has no sub-features\n\n"
                                )
                                numWarning += 1

                        elif featLvl2.type in [
                            "exon",
                            "CDS",
                            "intron",
                            "Shine_Dalgarno_sequence",
                        ]:
                            if (
                                featLvl2.location.start < featLvl1.location.start
                                or featLvl2.location.start > featLvl1.location.end
                            ):
                                errorMessage += (
                                    "Error: %s start falls outside the boundary of parent gene\n  %s gene Start: %d  -- End: %d\n  %s %s Start: %d\n\n"
                                    % (
                                        featLvl2.type,
                                        featLvl1.id,
                                        featLvl1.location.start,
                                        featLvl1.location.end,
                                        featLvl2.id,
                                        featLvl2.type,
                                        featLvl2.location.start,
                                    )
                                )
                                numError += 1
                            elif (
                                featLvl2.location.end < featLvl1.location.start
                                or featLvl2.location.end > featLvl1.location.end
                            ):
                                errorMessage += (
                                    "Error: %s end falls outside the boundary of parent gene\n  %s gene Start: %d  -- End: %d\n  %s %s End: %d\n\n"
                                    % (
                                        featLvl2.type,
                                        featLvl1.id,
                                        featLvl1.location.start,
                                        featLvl1.location.end,
                                        featLvl2.id,
                                        featLvl2.type,
                                        featLvl2.location.end,
                                    )
                                )
                                numError += 1
                            elif (
                                featLvl2.location.start != featLvl1.location.start
                                and featLvl2.location.end != featLvl1.location.end
                            ):
                                errorMessage += (
                                    "Error: %s does not match expected start/end pair based on parent gene boundary\n  %s gene Start: %d  -- End: %d\n  %s %s Start: %d  -- End: %d\n"
                                    % (
                                        featLvl2.type,
                                        featLvl1.id,
                                        featLvl1.location.start,
                                        featLvl1.location.end,
                                        featLvl2.id,
                                        featLvl2.type,
                                        featLvl2.location.start,
                                        featLvl2.location.end,
                                    )
                                )
                                numError += 1
                            errorMessage = (
                                "Error: %s sub-feature paired to gene and not mRNA feature\n\n"
                                % featLvl2.type
                            ) + errorMessage
                            numError += 1
                    if foundMRNA > 1:
                        errorMessage = (
                            "Error: Gene %s has more than one mRNA sub-feature\n\n"
                            % (featLvl1.id)
                        ) + errorMessage
                        numError += 1
                    elif foundMRNA == 0:
                        errorMessage = (
                            "Warning: Gene %s has no mRNA Feature\n\n" % (featLvl1.id)
                        ) + errorMessage
                        numWarning += 1

        locations.sort(key=lambda x: x.location.start)

        i = 0
        allowOver = 60

        while locations and i < len(locations) - 1:
            j = i + 1
            # To-Do: Check if i/j is tRNA and change below to allow flexibility
            if j + 1 < len(locations):
                if locations[i].type == "tRNA" or locations[j].type == "tRNA":
                    allowOver = 10
                else:
                    errorMessage = (
                        "Warning: %s mRNA has no CDS features\n\n" % (featLvl2.id)
                    ) + errorMessage
                    numWarning += 1

                if exons > 2:
                    errorMessage = (
                        "Error: %s mRNA has %d exons, expected no more than 2\n\n"
                        % (featLvl2.id, exons)
                    ) + errorMessage
                    numError += 1
                elif exons > 1 and len(exonList) > 1:
                    if not (
                        exonList[0].location.start > exonList[1].location.end
                        or exonList[1].location.start > exonList[0].location.end
                    ):
                        errorMessage = (
                            "Error: mRNA has 2 exons, and they don't match the expected start/end boundaries\n  mRNA %s Start: %d -- End: %d\n  Exon %s Start: %d -- End: %d\n  Exon %s Start: %d -- End: %d\n  Expected one exon's start to match the mRNA's start, the other to match its end, and neither to overlap.\n\n"
                            % (
                                featLvl2.id,
                                featLvl2.location.start,
                                featLvl2.location.end,
                                exonList[0].id,
                                exonList[0].location.start,
                                exonList[0].location.end,
                                exonList[1].id,
                                exonList[1].location.start,
                                exonList[1].location.end,
                            )
                        ) + errorMessage
                        numError += 1
                elif exons == 0:
                    errorMessage = (
                        "Warning: %s mRNA has no exons\n\n" % (featLvl2.id)
                    ) + errorMessage
                    numWarning += 1
                #            if(shines > 0):
                #              errorMessage = ("Error: %s mRNA has %d Shine Dalgarno sequences, expected no more than 1\n\n" % (featLvl2.id, shines))  + errorMessage

                if exons + CDS + shines == 0:
                    errorMessage = "Warning: %s mRNA has no sub-features\n\n"
                    numWarning += 1

            elif featLvl2.type in ["exon", "CDS", "intron", "Shine_Dalgarno_sequence"]:
                if (
                    featLvl2.location.start < featLvl1.location.start
                    or featLvl2.location.start > featLvl1.location.end
                ):
                    errorMessage += (
                        "Error: %s start falls outside the boundary of parent gene\n  %s gene Start: %d  -- End: %d\n  %s %s Start: %d\n\n"
                        % (
                            featLvl2.type,
                            featLvl1.id,
                            featLvl1.location.start,
                            featLvl1.location.end,
                            featLvl2.id,
                            featLvl2.type,
                            featLvl2.location.start,
                        )
                    )
                    numError += 1
                elif (
                    featLvl2.location.end < featLvl1.location.start
                    or featLvl2.location.end > featLvl1.location.end
                ):
                    errorMessage += (
                        "Error: %s end falls outside the boundary of parent gene\n  %s gene Start: %d  -- End: %d\n  %s %s End: %d\n\n"
                        % (
                            featLvl2.type,
                            featLvl1.id,
                            featLvl1.location.start,
                            featLvl1.location.end,
                            featLvl2.id,
                            featLvl2.type,
                            featLvl2.location.end,
                        )
                    )
                    numError += 1
                elif (
                    featLvl2.location.start != featLvl1.location.start
                    and featLvl2.location.end != featLvl1.location.end
                ):
                    errorMessage += (
                        "Error: %s does not match expected start/end pair based on parent gene boundary\n  %s gene Start: %d  -- End: %d\n  %s %s Start: %d  -- End: %d\n"
                        % (
                            featLvl2.type,
                            featLvl1.id,
                            featLvl1.location.start,
                            featLvl1.location.end,
                            featLvl2.id,
                            featLvl2.type,
                            featLvl2.location.start,
                            featLvl2.location.end,
                        )
                    )
                    numError += 1
                errorMessage = (
                    "Error: %s sub-feature paired to gene and not mRNA feature\n\n"
                    % featLvl2.type
                ) + errorMessage
                numError += 1
            if foundMRNA > 1:
                errorMessage = (
                    "Error: Gene %s has more than one mRNA sub-feature\n\n"
                    % (featLvl1.id)
                ) + errorMessage
                numError += 1
            elif foundMRNA == 0:
                errorMessage = (
                    "Warning: Gene %s has no mRNA Feature\n\n" % (featLvl1.id)
                ) + errorMessage
                numWarning += 1

    locations.sort(key=lambda x: x.location.start)

    i = 0
    allowOver = 60
    while locations and i < len(locations) - 1:
        j = i + 1
        # To-Do: Check if i/j is tRNA and change below to allow flexibility
        if j + 1 < len(locations):
            if locations[i].type == "tRNA" or locations[j].type == "tRNA":
                allowOver = 10
            else:
                allowOver = 60

        while j < len(locations):
            if (locations[i].location.end - locations[j].location.start) >= allowOver:
                errorMessage = (
                    "Error: %s %s overlaps %s %s by %d or more bases\n  %s %s Start: %d -- End: %d\n  %s %s Start: %d -- End: %d\n\n"
                    % (
                        locations[i].type,
                        locations[i].id,
                        locations[j].type,
                        locations[j].id,
                        allowOver,
                        locations[i].type,
                        locations[i].id,
                        locations[i].location.start,
                        locations[i].location.end,
                        locations[j].type,
                        locations[j].id,
                        locations[j].location.start,
                        locations[j].location.end,
                    )
                ) + errorMessage
                numError += 1
            elif (
                locations[i].location.end - locations[j].location.start
            ) > 0 and allowOver == 10:
                errorMessage = (
                    "Warning: %s %s overlaps %s %s, may wish to verify\n  %s %s Start: %d -- End: %d\n  %s %s Start: %d -- End: %d\n\n"
                    % (
                        locations[i].type,
                        locations[i].id,
                        locations[j].type,
                        locations[j].id,
                        locations[i].type,
                        locations[i].id,
                        locations[i].location.start,
                        locations[i].location.end,
                        locations[j].type,
                        locations[j].id,
                        locations[j].location.start,
                        locations[j].location.end,
                    )
                )
                numWarning += 1
            j += 1

        i += 1

    if errorMessage:
        out_errorlog.write(errorMessage)

    out_errorlog.write(
        "Finished with %d Errors and %d Warnings" % (numError, numWarning)
    )

    # For each terminator/tRNA
    # Bother Cory later


#  print(dir(record.features))
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Validate that contents of GFF3 are correct"
    )
    parser.add_argument("gff3In", type=argparse.FileType("r"), help="GFF3 source file")
    parser.add_argument(
        "--out_errorlog",
        type=argparse.FileType("w"),
        help="Output Error Log",
        default="test-data/errorlog.txt",
    )
    args = parser.parse_args()
    table_annotations(**vars(args))
