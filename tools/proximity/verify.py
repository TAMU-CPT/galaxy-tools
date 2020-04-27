from Bio import SeqIO
import re
from Bio.Blast import NCBIXML

def readGBK(files,search_list):
    if files:
        for idx, file in enumerate(files):

            if idx == 0:
                print("Parsing - "+file)
                record = SeqIO.read(file, "genbank")
                gbk_matches = []
                for feature in record.features:
                    try:
                        #print("did find it =_=_="+str(feature.type))
                        if searchInput(str(feature.qualifiers["product"]),search_list=search_list) or searchInput(str(feature.qualifiers["note"]),search_list=search_list) or search_list(str(feature.qualifiers["dbxref"])):
                            print(feature.qualifiers["dbxref"])
                            print(feature.qualifiers["product"])
                            print(feature.qualifiers["note"])
                            gbk_matches.extend(searchInput(str(feature),search_list=search_list))
                    except KeyError:
                        #print("didn't find it _=_=_="+str(feature.type))
                        #print(feature.type)
                        continue
                gbk_matches = list(set(gbk_matches))
                print(gbk_matches)
            else:
                print("Parsing - "+file)
                record = SeqIO.read(file, "genbank")
                for feature in record.features:
                    print(feature)
                    gbk_matches.extend(searchInput(str(feature),search_list=search_list))
                gbk_matches = list(set(gbk_matches))

        return gbk_matches
    else:
        pass

def searchInput(input, search_list,blast=False,q_id=""):
    """ Takes an input search string, and returns uniques of passing """
    output = []
    for search_term in search_list:
        if blast:
            #print(input)
            if re.search(re.escape(search_term), input, flags=re.IGNORECASE):
                add_query = "QueryID: "+str(q_id)+"\nSearchQuery: "+search_term+"\nMatch: "+input+"\n"
                output.extend([add_query])
            else:
                continue
        #print(search_term)
        #st = r"\b"+search_term+r"\b"
        else:
            if re.search(re.escape(search_term), input,flags=re.IGNORECASE):
                #print(search_term+" -> was found")
                output.extend([input])
            else:
                continue
    return list(set(output))

def readBLAST(files,search_list):
    if files:
        for idx, file in enumerate(files):
            if idx == 0:
                print("Parsing - "+file)
                blast_records = NCBIXML.parse(open(file))
                blast_matches = []
                for blast_record in blast_records:
                    #print("********")
                    #print(blast_record.query)
                    for desc in blast_record.descriptions:
                        #print("====")
                        #print(desc)
                        #print(feature.id)
                        pretty = prettifyXML(str(desc))
                        for each_ret in pretty:
                            blast_matches.extend(searchInput(each_ret,search_list=search_list,blast=True,q_id=blast_record.query))
                blast_matches = list(set(blast_matches))
                print(blast_matches)
            else:
                print("Parsing - "+file)
                blast_records = NCBIXML.parse(open(file))
                for blast_record in blast_records:
                    for desc in blast_record.descriptions:
                        pretty = prettifyXML(str(desc))
                        blast_matches.extend(searchInput(each_ret,search_list=search_list,blast=True))
                blast_matches = list(set(blast_matches))
            return blast_matches
    else:
        pass

######## prettify-XML function
def prettifyXML(input):
    """ prettifies a string input from a BLAST-xml """
    s = input
    split = s.split(">")

    return split



if __name__ == "__main__":
    print("import worked")
    #readGBK(["test-data/NC_008720-N4-Apollo.genbank"],["endolysin","holin"])
    b = readBLAST(["test-data/Galaxy27-NR.blastxml.xml"],["pyrophosphohydrolase"])
    for i in b:
        print(i)