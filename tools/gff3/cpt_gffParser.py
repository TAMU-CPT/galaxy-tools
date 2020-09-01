from Bio.SeqFeature import FeatureLocation, CompoundLocation
from Bio import SeqIO, SeqFeature, SeqRecord
from Bio.Seq import Seq

disallowArray = ["&", ",", ";", "="]
validArray = ["%26", "%2C", "%3B", "%3D"]

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
      return "", None
    if line[0] == "#":
      if len(line) > 2 and line[1] == "#":
        return None, line[2:-1], None 
      # else handle ## Pragmas
      else: 
        return None, None, None
        

    errorMessage = ""

    fields = line.split("\t")
    if len(fields[0]) != 9:
      errorMessage += "GFF3 is a 9-column tab-separated format, line has " + str(len(fields)) + " columns.\n"
      if len(fields[0]) > 9:
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
      

    
    if fields[7] not in ['.', '0', '1', '2']:
      errorMessage += "Expected 0, 1, 2, or . for Phase field value, actual value is '" + fields[7] + "'.\n"
    elif fields[7] =='.' and fields[1] == "CDS":
      errorMessage += "Expected 0, 1, or 2 in Phase field for CDS-type feature, actual value is '" + fields[7] + "'.\n"

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
        if len(qualDict[x]) > 1:
          errorMessage += "More than one ID supplied for feature.\n"
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
        
def gffParse(gff3In):
    fastaDirective = False
    errOut = ""
    featList = []
    featInd = 0
    seekingParent = []
    lineInd = 0
    orgDict = {}
    seekParentDict = {}
    indDict = {}

    for line in gff3In:
      lineInd += 1
      err, prag, res = lineAnalysis(line)
      if err:
        errOut += (str(lineInd) + ": " + err + "\n")
      if prag and not res:
        if prag == "FASTA":
          fastaDirective = True
      if res:
        if prag not in indDict.keys():
          indDict[prag] = 0
          orgDict[prag] = []
          seekParentDict[prag] = []
        else:
          indDict[prag] += 1
        orgDict[prag].append(res)
        if "Parent" in res.qualifiers.keys():
          seekParentDict[prag].append(indDict[prag])
    if errOut:
      raise Exception("Failed GFF Feature Parsing with: \n" + errOut)
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
      for i in orgDict[x]:
        if "Parent" not in i.qualifiers.keys():
          finalOrgHeirarchy.append(i)
      res.append(SeqRecord.SeqRecord(None, x, "<unknown name>", "<unknown description>", None, finalOrgHeirarchy, None, None))
    return res
