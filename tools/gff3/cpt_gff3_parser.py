#!/usr/bin/env python

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import argparse

disallowArray = ["&", ",", ";", "="]
validArray = ["%26", "%2C", "%3B", "%3D"]

validID = '.:^*$@!+_?-|'

def validateID(idIn):
    badChar = []
    for x in idIn:
      if (ord(x) > 47 and ord(x) < 58) or (ord(x) > 64 and ord(x) < 91) or (ord(x) > 96 and ord(x) < 123) or (x in validID):
        continue
      else:
        if not(x in badChar):
          badChar.append(x)
    return badChar 

def replaceBadChars(qualIn):
    newQual = ""
    for x in qualIn:
      goodVal = True
      for y in range(0, len(disallowArray)):
        if x == disallowArray[y]:
          goodVal = False
          newQual += validArray[y]
      if goodVal:
        newQual += x
    return newQual

def validateQual(qualIn):
    badChar = []
    for x in qualIn:
      if x in disallowArray:
        if not(x in badChar):
          badChar.append(x)
    return badChar 

def encodeFromLookahead(remLine):
    for x in remLine:
      if x == "=":
        return False
      if x == ";" or x == ",":
        return True
    return True

def lineAnalysis(line):
    IDName = ""
    startLoc = -1
    endLoc = -1
    score = None
    if len(line) == 0:
      return "", None
    if line[0] == "#":
      if len(line) > 1 and line[1] != "#":
        return "", None 
      # else handle ## Pragmas
      else:
        return "", None
    errorMessage = ""

    fields = line.split("\t")
    if len(fields) != 9:
      errorMessage += "GFF3 is a 9-column tab-separated format, line has " + str(len(fields)) + " columns.\n"
      if len(fields) > 9:
        errorMessage += "Possible unescaped tab in a qualifier field.\n"
        return errorMessage, None

    for x in range(0, len(fields)):
      if fields[x] == "":
        errorMessage += "Field #" + str(x + 1) + " is empty. Please supply correct or default value.\n"
    if errorMessage != "":
      return errorMessage, None

    idEval = validateID(fields[0])
    if len(idEval) != 0:
      errorMessage += "Organism ID contains the following invalid characters: " + str(idEval) + ".\n"

    # fields[1]
    
    # fields[2]

    isNum = True
    for x in fields[3]:  
      if not(ord(x) > 47 and ord(x) < 58):
        errorMessage += "Feature location start contains non-numeric character.\n"
        isNum = False
        break
    if isNum:
      startLoc = int(fields[3])
    
    isNum = True
    for x in fields[4]:
      if not(ord(x) > 47 and ord(x) < 58):
        errorMessage += "Feature location end contains non-numeric character.\n"
        isNum = False
        break
    if isNum:
      endLoc = int(fields[4])

    if endLoc < startLoc:
      errorMessage += "Feature Location end is less than start (GFF  spec requires all features, regardless of strand, to have the lower number as the start).\n"

    # fields[5]
    isNum = False
    foundDot = False
    if fields[5] != ".":
      isNum = True
      for x in fields[5]:
        if not(ord(x) > 47 and ord(x) < 58):
          if x == "." and not foundDot:
            foundDot = True
          else:
            errorMessage += "Feature score is a non-numeric structure.\n"
            isNum = False
            break
      if isNum:
        score = float(fields[5])


    if len(fields[6]) != 1 or (not(fields[6] in '-+.?')):
      errorMessage += "Feature strand must be '+', '-', '.', or '?', actual value is '" + fields[6] + "'.\n"
      


    if fields[7] != '.' and fields[2] != "CDS":
      errorMessage += "Expected '.' in Phase field for non-CDS feature, actual value is '" + fields[7] + "'.\n"
    elif fields[2] == "CDS":
      for x in fields[7]:
        if not(ord(x) > 47 and ord(x) < 58):
          errorMessage += "Non-numeric value in Phase field for CDS type feature '" + fields[7] + "'.\n"

    keyName = ""
    valNames = [""]
    valInd = 0
    parseMode = 0  
    qualDict = {} 
    for x in range(0, len(fields[8])):
      currChar = fields[8][x]
      if parseMode == 0:
        if not (currChar in "=,;"):
          keyName += currChar
        elif currChar == "=":
          if len(keyName) == 0:
            errorMessage += "No ID name supplied for a value in the qualifiers field, aborting.\n"
            break
          parseMode = 1
          continue
        else: #Encode special char
          keyName += "%" + str(hex(ord(currChar)))
      elif parseMode == 1:
        if not (currChar in "=,;"):
          valNames[valInd] += currChar
        elif currChar == ",":
          valInd += 1
          valNames.append("")
        elif currChar == "=":
          valNames[valInd] += "%3D"
        else:
          if x == len(fields[8]) - 1: # Assume if last char in fields[8] is a semicolon, then just the end of qualifier 
            parseMode = 2
          elif encodeFromLookahead(fields[8][x+1:]):
            valNames[valInd] += "%3B"
            continue
          else:
            parseMode = 2
      if parseMode == 2: # Do not elif, this is used as a wrapup for each qualifier and we want it checked if parsemode == 1 incremented itself 
        if not qualDict.keys() or keyName not in list(qualDict.keys()):
          qualDict[keyName] = valNames
        else:
          for x in valNames:
            qualDict[keyName].append(x)
        keyName = ""      
        valNames = [""]
        valInd = 0
        parseMode = 0

    for x in list(qualDict.keys()):
      if x == "ID":
        if len(qualDict[x]) > 1:
          errorMessage += "More than one ID supplied for feature.\n"
        IDName = qualDict[x][0]

    if startLoc == -1 or endLoc == -1 or (not(fields[6] in '-+.?')):
      errorMessage += "Unable to construct feature location, aborting.\n"
    elif fields[6] == '+':
      featLoc = FeatureLocation(startLoc - 1, endLoc - 1, strand = +1)
    elif fields[6] == '-':
      featLoc = FeatureLocation(startLoc - 1, endLoc - 1, strand = -1)
    else:
      featLoc = FeatureLocation(startLoc - 1, endLoc - 1, strand = 0)
     
    

    if errorMessage != "":
      return errorMessage, None
      
    return None, SeqFeature(featLoc, fields[2], '', featLoc.strand, IDName, qualDict, None, None, None)   

def parseGFF(fileIn):
    lineIn = fileIn.readline()
    lineNum = 1
    lookupParents = []
    loadedFeats = []
    chainFeats = []
    chainFeats.append([])
    subFeatIndices = []
    finalizeList = []
    foundIDs = []
    highestLoc = -1
    while lineIn:
      errOut, featOut = lineAnalysis(lineIn)
      if errOut:
        print(str(lineNum) + ": " + errOut)
      elif featOut:
        if "ID" in list(featOut.qualifiers.keys()): 
          if featOut.qualifiers["ID"][0] in foundIDs:
            print("Duplicate IDs: " + featOut.qualifiers["ID"][0] + ", aborting.\n")
            exit()
          else:
            foundIDs.append(featOut.qualifiers["ID"][0])

        if "Parent" in list(featOut.qualifiers.keys()):
          lookupParents.append(len(loadedFeats))
          finalizeList.append(-1)
        else: 
          finalizeList.append(0)
        loadedFeats.append(featOut)
        subFeatIndices.append([])
        

        highestLoc = int(max(highestLoc, max(featOut.location.start, featOut.location.end)))
      lineIn = fileIn.readline()
      lineNum += 1
    depth = 0
    for child in lookupParents:
      foundParent = False
      lookupID = loadedFeats[child].qualifiers["Parent"][0]
      for i in range(0, len(loadedFeats)):
        if loadedFeats[i].id == lookupID:
          foundParent = True
          if finalizeList[i] < 0:
            lookupParents.append(child)
          else:
            finalizeList[child] = finalizeList[i] + 1
            subFeatIndices[child].append(i)
            depth = int(max(depth, finalizeList[child]))
      if not foundParent:
        print("Parent ID " + lookupID + " of sub-feature " + loadedFeats[child].id + " not found, aborting.")
        exit()
      
    depthI = 0
    while depthI < depth: 
      chainFeats.append([])
      depthI += 1

    for x in range(0, len(finalizeList)):
      chainFeats[finalizeList[x]].append(x)
    

    depthI = depth - 1
    featNum = 0
    while depthI >= 0:
      for x in chainFeats[depthI]:
        subFeats = []
        for y in subFeatIndices[x]:
          subFeats.append(loadedFeats[y])
        #print("Make Feat")
        loadedFeats[x] = SeqFeature(loadedFeats[x].location, loadedFeats[x].type, '', loadedFeats[x].strand, loadedFeats[x].id, loadedFeats[x].qualifiers, subFeats, None, None)
      depthI -= 1
    
    outRec = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF", IUPAC.protein), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=None, features=loadedFeats, annotations=None, letter_annotations=None)   

    print(outRec)
          
          
      
    #print(loadedFeats)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Validate that contents of GFF3 are correct"
    )
    parser.add_argument("gff3In", type=argparse.FileType("r"), help="GFF3 source file")
    args = parser.parse_args()   
    parseGFF(args.gff3In)
         
