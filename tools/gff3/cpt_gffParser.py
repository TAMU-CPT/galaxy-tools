from Bio.SeqFeature import FeatureLocation, CompoundLocation
from Bio import SeqIO, SeqFeature, SeqRecord
from Bio.Seq import Seq, UnknownSeq
import sys
import urllib

#import urlparse.urlparse
#import urllib.parse

disallowArray = ["&", ",", ";", "="]
validArray = ["%26", "%2C", "%3B", "%3D"]
encoders = "ABCDEFGHIJKLMNOPQRSTUVWXYZ01234567890"

validID = '.:^*$@!+_?-|'

class gffSeqFeature(SeqFeature.SeqFeature):
    def __init__(
        self,
        location=None,
        type="",
        location_operator="",
        strand=None,
        id="<unknown id>",
        qualifiers=None,
        sub_features=None,
        ref=None,
        ref_db=None,
    ):
        """Reimplementation of SeqFeature for use with GFF3 Parsing
        Does not remove the sub_feature functionality, as unlike
        Genbank, this is baked into the core concept of GFF
        """
        if (
            location is not None
            and not isinstance(location, FeatureLocation)
            and not isinstance(location, CompoundLocation)
        ):
            raise TypeError(
                "FeatureLocation, CompoundLocation (or None) required for the location"
            )
        self.location = location
        self.type = type
        if location_operator:
            # TODO - Deprecation warning
            self.location_operator = location_operator
        if strand is not None:
            # TODO - Deprecation warning
            self.strand = strand
        self.id = id
        if qualifiers is None:
            qualifiers = OrderedDict()
        self.qualifiers = qualifiers
        if sub_features is None:
            sub_features = []
        self._sub_features = sub_features
        if ref is not None:
            # TODO - Deprecation warning
            self.ref = ref
        if ref_db is not None:
            # TODO - Deprecation warning
            self.ref_db = ref_db

    def _get_subfeatures(self):
        """Get function for the sub_features property (PRIVATE)."""
        try:
            return self._sub_features
        except AttributeError:
            return None

    def _set_subfeatures(self, value):
        """Set function for the sub_features property (PRIVATE)."""
        if value:
            if isinstance(value, list):
                self._sub_features = value
            else:
                raise ValueError("sub_feature must be a list of gffSeqFeature objects")

    sub_features = property(
        fget=_get_subfeatures,
        fset=_set_subfeatures,
        doc="Sub-features for GFF Heirarchy",
    )


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
    if len(line) == 0 or line == "\n":
      return None, None, None
    if line[0] == "#":
      if len(line) > 2 and line[1] == "#":
        return None, line[2:-1], None 
      # else handle ## Pragmas
      else: 
        return None, None, None
        

    errorMessage = ""

    fields = line.split("\t")
    if len(fields) != 9:
      errorMessage += "GFF3 is a 9-column tab-separated format, line has " + str(len(fields)) + " columns.\n"
      if len(fields) > 9:
        errorMessage += "Possible unescaped tab in a qualifier field.\n"
        return errorMessage, None, None

    for x in range(0, len(fields)):
      if fields[x] == "":
        errorMessage += "Field #" + str(x + 1) + " is empty. Please supply correct or default value.\n"
    if errorMessage != "":
      return errorMessage, None, None

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
    if fields[5] != ".":
      try:
        score = float(fields[5])
      except:
        score = None
        errorMessage += "Score field could not be interpreted as a floating-point (real) number. Ensure notation is correct.\n"

    if len(fields[6]) != 1 or (not(fields[6] in '-+.?')):
      errorMessage += "Feature strand must be '+', '-', '.', or '?', actual value is '" + fields[6] + "'.\n"
      

    
    if fields[7] not in ['.', '0', '1', '2']:
      errorMessage += "Expected 0, 1, 2, or . for Phase field value, actual value is '" + fields[7] + "'.\n"
    elif fields[7] =='.' and fields[1] == "CDS":
      errorMessage += "Expected 0, 1, or 2 in Phase field for CDS-type feature, actual value is '" + fields[7] + "'.\n"

    keyName = ""
    valNames = [""]
    valInd = 0
    parseMode = 0  
    qualDict = {} 
    contCounter = 0
    for x in range(0, len(fields[8])):
      currChar = fields[8][x]
      if contCounter:
        contCounter += -1
        continue
      if parseMode == 0:
        if not (currChar in "=,;%"):
          keyName += currChar
        elif currChar == "=":
          if len(keyName) == 0:
            errorMessage += "No ID name supplied for a value in the qualifiers field, aborting.\n"
            break
          parseMode = 1
          continue
        elif currChar == "%" and (fields[8][x+1] in encoders) and (fields[8][x+2] in encoders):
          keyName += urllib.unquote(fields[8][x:x+3])
          contCounter = 2
        else: #Encode special char
          encoded = str(hex(ord(currChar)))
          keyName += "%" + encoded[2:].upper()
      elif parseMode == 1:
        if not (currChar in "=,;%\n"):
          valNames[valInd] += currChar
        elif currChar == ",":
          valInd += 1
          valNames.append("")
        elif currChar == "=":
          valNames[valInd] += "%3D"
        elif currChar == "%" and (fields[8][x+1] in encoders) and (fields[8][x+2] in encoders):
          valNames[valInd] += urllib.unquote(fields[8][x:x+3])
          contCounter = 2
        else:
          if x == len(fields[8]) - 1: # Assume if last char in fields[8] is a semicolon, then just the end of qualifier 
            parseMode = 2
          elif encodeFromLookahead(fields[8][x+1:]):
            valNames[valInd] += "%3B"
            continue
          else:
            parseMode = 2
      if parseMode == 2: # Do not elif, this is used as a wrapup for each qualifier and we want it checked if parsemode == 1 incremented itself 
        if keyName not in qualDict.keys():
          qualDict[keyName] = valNames
        else:
          for x in valNames:
            qualDict[keyName].append(x)
        keyName = ""      
        valNames = [""]
        valInd = 0
        parseMode = 0

    for x in qualDict.keys():
      if x == "ID":
        #if len(qualDict[x]) > 1:
          #errorMessage += "More than one ID supplied for feature.\n"
        IDName = qualDict[x][0]

    if startLoc == -1 or endLoc == -1 or (not(fields[6] in '-+.?')):
      errorMessage += "Unable to construct feature location, aborting.\n"
    elif fields[6] == '+':
      featLoc = FeatureLocation(startLoc - 1, endLoc, strand = +1)
    elif fields[6] == '-':
      featLoc = FeatureLocation(startLoc - 1, endLoc, strand = -1)
    else:
      featLoc = FeatureLocation(startLoc - 1, endLoc, strand = 0)
     
    

    if errorMessage != "":
      return errorMessage, None, None
      
    return None, fields[0], gffSeqFeature(featLoc, fields[2], '', featLoc.strand, IDName, qualDict, None, None, None)   
        
def gffParse(gff3In, base_dict = {}, outStream = sys.stderr):
    fastaDirective = False
    errOut = ""
    featList = []
    featInd = 0
    seekingParent = []
    lineInd = 0
    orgDict = {}
    seekParentDict = {}
    indDict = {}
    seqDict = {}
    regionDict = {}
    currFastaKey = ""
    

    for line in gff3In:
      lineInd += 1
      err = None
      prag = None
      res = None
      if line[0] == ">":	# For compatibility with Artemis-style GFF
        fastaDirective = True
      if not fastaDirective:
        err, prag, res = lineAnalysis(line)
      else:
        if line[0] == ">":
          currFastaKey = line[1:-1]
        elif line[0] == "#":
          continue
        elif line:
          seqDict[currFastaKey] += (line[:-1]).strip()
      if err:
        errOut += ("Line " + str(lineInd) + ": " + err + "\n")
      if prag and not res:
        prag = prag.split(" ")
        if prag[0] == "FASTA":
          fastaDirective = True
        elif prag[0] == "sequence-region":
          regionDict[prag[1]] = [int(prag[2]) - 1, int(prag[3])]
      if res:
        if prag not in indDict.keys():
          indDict[prag] = 0
          orgDict[prag] = []
          seekParentDict[prag] = []
          seqDict[prag] = ""
        else:
          indDict[prag] += 1
        orgDict[prag].append(res)
        if "Parent" in res.qualifiers.keys():
          seekParentDict[prag].append(indDict[prag])
    if errOut:
      outStream.write(errOut + "\n")
      raise Exception("Failed GFF Feature Parsing, error log output to stderr")
    for org in seekParentDict.keys():
      for ind in seekParentDict[org]:
        for x in orgDict[org][ind].qualifiers['Parent']:
          for y in orgDict[org]:
            found = False
            if y.id == x:
              y.sub_features.append(orgDict[org][ind])
              found = True
              break
          if not found:
            raise Exception("Unable to find parent " + x + " of feature " + orgDict[org][ind].id)

    res = []
    for x in orgDict.keys():
      finalOrgHeirarchy = []
      annoteDict = {}
      for i in orgDict[x]:
        if "Parent" not in i.qualifiers.keys():
          finalOrgHeirarchy.append(i)
      if seqDict[x]:
        if x in regionDict.keys():
          annoteDict["sequence-region"] = (x, regionDict[x][0], regionDict[x][1])
          if len(seqDict[x]) < regionDict[x][1] - regionDict[x][0]:
            seqDict[x] += "?" * (regionDict[x][1] - regionDict[x][0] - len(seqDict[x]))
          else:
            seqDict[x] = seqDict[x][regionDict[x][0]:regionDict[x][1]]
        else:
          annoteDict["sequence-region"] = (x, 0, len(seqDict[x]))
        seqDict[x] = Seq(seqDict[x])
      elif x in regionDict.keys():
        annoteDict["sequence-region"] = (x, regionDict[x][0], regionDict[x][1])
        seqDict[x] = UnknownSeq(regionDict[x][1] - regionDict[x][0])
      else:
        seqDict[x] = None
      res.append(SeqRecord.SeqRecord(seqDict[x], x, "<unknown name>", "<unknown description>", None, finalOrgHeirarchy, annoteDict, None))
    standaloneList = []
    for x in regionDict.keys():
      if x not in orgDict.keys():
        standaloneList.append(x)
    """for x in standaloneList:
      annoteDict = {}
      annoteDict["sequence-region"] = (x, regionDict[x][0], regionDict[x][1])
      if x in seqDict.keys():
        if len(seqDict[x]) < regionDict[x][1] - regionDict[x][0]:
          seqDict[x] += "?" * (regionDict[x][1] - regionDict[x][0] - len(seqDict[x]))
        else:
          seqDict[x] = seqDict[x][regionDict[x][0]:regionDict[x][1]]
      else:
        seqDict[x] = UnknownSeq(regionDict[x][1] - regionDict[x][0])
      res.append(SeqRecord.SeqRecord(seqDict[x], x, "<unknown name>", "<unknown description>", None, None, annoteDict, None))
    """
    for x in base_dict.keys():
      found = False
      for y in res:
        if x == y.id:
          found = True
          y.name = base_dict[x].name
          y.description = base_dict[x].description
          y.seq = base_dict[x].seq
          break
      

    return res

def printFeatLine(inFeat, orgName, source = 'feature', score = None, shift = None, outStream = sys.stdout):
    line = orgName + "\t"
    if source:
      line += source + "\t"
    else:
      line += ".\t" 
    line += inFeat.type + "\t"
    line += str(min(inFeat.location.start, inFeat.location.end) + 1) + "\t" + str(max(inFeat.location.start, inFeat.location.end)) + "\t"
    if score:
      line += str(score) + "\t"
    else:
      line += ".\t"
    if inFeat.location.strand == None:
      line += ".\t"
    elif inFeat.location.strand == 1:
      line += "+\t"
    elif inFeat.location.strand == -1:
      line += "-\t"
    else: 
      line += "?\t"
    if shift:
      line += str(shift) + "\t"
    elif inFeat.type == "CDS":
      line += "0\t"
    else:
      line += ".\t"
    for qual in inFeat.qualifiers.keys():
      for keyChar in str(qual):
        if keyChar in "%,=;":
          encoded = str(hex(ord(keyChar)))
          line += "%" + encoded[2:].upper()
        else:
          line += keyChar
      line += "="
      for ind in range(0, len(inFeat.qualifiers[qual])):
        for valChar in str(inFeat.qualifiers[qual][ind]):
          if valChar in "%,=;":
            encoded = str(hex(ord(valChar)))
            line += "%" + encoded[2:].upper()
          else:
            line += valChar
        if ind < len(inFeat.qualifiers[qual]) - 1:
          line += ","
        else:
          line += ";"
    outStream.write(line + "\n")
  
    if type(inFeat) == gffSeqFeature and inFeat.sub_features: 
      for x in inFeat.sub_features:
        printFeatLine(x, orgName, source, score, shift, outStream)

def gffWrite(inRec, outStream = None):
    if not outStream:
      outStream = sys.stdout
    outStream.write("##gff-version 3\n")
    if not inRec:
      return
    if type(inRec) != list:
      inRec = [inRec]
    for rec in inRec:
      if "sequence-region" in rec.annotations.keys():
        outStream.write("##sequence-region " + rec.annotations["sequence-region"][0] + " " + str(rec.annotations["sequence-region"][1] + 1) + " " + str(rec.annotations["sequence-region"][2]) + "\n")
      #else make one up based on feature locations?
      elif rec.seq:
        outStream.write("##sequence-region " + rec.id + " 1 " + str(len(rec.seq)) +"\n")
      for feat in rec.features:
          printFeatLine(feat, rec.id, outStream = outStream)        
