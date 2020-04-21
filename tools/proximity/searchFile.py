##### User input File(s), that are BLAST XML, gff3, and/or Genbank and then searched for containing user designated terms

import argparse 
import explodeJSON as ej
import gffutils # THIS IS REQUIREMENT
from Bio.Blast import NCBIXML
from Bio import SeqIO
import re
import os

####### TERM FUNCTIONS
def dbaseTerms(terms,galaxy=True):
    """ Index into dictionary object and retrieve all desired terms """
    if galaxy:
        db_path = "/galaxy/tools/cpt2/galaxy-tools/tools/proximity/data/lysis-family-expanded.json"
    else:
        db_path = "data/lysis-family-expanded.json"
    db = ej.explodeJSON(db_path)
    db = db.readJSON()
    dbase_terms = []
    if terms:
        for term in terms:
            index_list = term.split(",")
            for t in index_list:
                dbase_terms.extend(db[t])
        return dbase_terms
    else:
        pass



def userTerms(file,text):
    """ Select terms input by user """
    user_terms = []
    if file:
        terms = open(file.name).read().splitlines()
        user_terms.extend(terms)
    else:
        pass
    if text:
        user_terms.extend(text)
    else:
        pass

    return user_terms


def glueTerms(dbase_terms, user_terms):
    """ glue dbaseTerms and userTerms together for eventual query item """
    glued = []
    if dbase_terms:
        glued.extend(dbase_terms)
    else:
        pass
    if user_terms:
        glued.extend(user_terms)
    else:
        pass

    return glued

####### FILE FUNCTIONS
def glueFiles(gff,gbk,fa,blast):
    """ glue giles into one list...I think this is a decent way to go about this...#CHECK LATER#... """
    files = []
    gffs = []
    gbks = []
    blasts = []
    if gff:
        gffs.extend(gff)
    else:
        pass
    if gbk:
        gbks.extend(gbk)
    else:
        pass
    fas = []
    if fa:
        fas.extend(fa)
    else:
        pass
    if blast:
        blasts.extend(blast)
    else:
        pass
    files = [gffs,gbks,fas,blasts]

    return files

######## PARSE FILE FUNCTIONS
def readGFF3(files,search_list):
    " Searches through gff3 file(s) and appends "
    if files:
        for idx, file in enumerate(files):
            if idx == 0:
                print("Parsing - "+file.name)
                db = gffutils.create_db(file.name,dbfn="file.db",force=True,keep_order=False)
                db = gffutils.FeatureDB("file.db")
                features = db.all_features()
                gff3_matches = []
                for feature in features:
                    gff3_matches.extend(searchInput(str(feature), search_list=search_list))
                gff3_matches = list(set(gff3_matches)) # make sure we don't fluff the list
            else:
                print("Parsing - "+file.name)
                db = gffutils.create_db(file.name,dbfn=str(idx)+"_file.db",force=True,keep_order=False)
                db = gffutils.FeatureDB(str(idx)+"_file.db")
                features = db.all_features()
                for feature in features:
                    gff3_matches.extend(searchInput(str(feature), search_list=search_list))
                gff3_matches = list(set(gff3_matches)) # make sure we don't fluff the list
        return gff3_matches
    else:
        pass

def readGBK(files,search_list):
    if files:
        for idx, file in enumerate(files):
            if idx == 0:
                print("Parsing - "+file.name)
                record = SeqIO.read(file.name, "genbank")
                gbk_matches = []
                for feature in record.features:
                    gbk_matches.extend(searchInput(str(feature),search_list=search_list))
                gbk_matches = list(set(gbk_matches))
            else:
                print("Parsing - "+file.name)
                record = SeqIO.read(file.name, "genbank")
                for feature in record.features:
                    gbk_matches.extend(searchInput(str(feature),search_list=search_list))
                gbk_matches = list(set(gbk_matches))
        return gbk_matches
    else:
        pass

def readFASTA(files,search_list):
    if files:
        for idx, file in enumerate(files):
            if idx == 0:
                print("Parsing - "+file.name)
                record = SeqIO.parse(file.name, "fasta")
                fa_matches = []
                for feature in record:
                    fa_matches.extend(searchInput(feature.description,search_list=search_list))
                fa_matches = list(set(fa_matches))
            else:
                print("Parsing - "+file.name)
                record = SeqIO.parse(file.name, "fasta")
                for feature in record:
                    fa_matches.extend(searchInput(feature.description,search_list=search_list))
                fa_matches = list(set(fa_matches))
        return fa_matches
    else:
        pass

def readBLAST(files,search_list):
    if files:
        for idx, file in enumerate(files):
            if idx == 0:
                print("Parsing - "+file.name)
                record = NCBIXML.parse(open(file.name))
                blast_matches = []
                for feature in record:
                    #print(feature.descriptions)
                    for desc in feature.descriptions:
                        blast_matches.extend(searchInput(str(desc),search_list=search_list))
                blast_matches = list(set(blast_matches))
            else:
                print("Parsing - "+file.name)
                record = NCBIXML.parse(open(file.name))
                for feature in record:
                    for desc in feature.descriptions:
                        blast_matches.extend(searchInput(str(desc),search_list=search_list))
                blast_matches = list(set(blast_matches))
            return blast_matches
    else:
        pass


######## SEARCH FILE FUNCTIONS
def searchInput(input, search_list):
    """ Takes an input search string, and returns uniques of passing """
    output = []
    for search_term in search_list:
        #print(search_term)
        #st = r"\b"+search_term+r"\b"
        if re.search(re.escape(search_term), input,flags=re.IGNORECASE):
            #print(search_term+" -> was found")
            output.extend([input])
        else:
            continue
    return list(set(output))


########## Output File Writer
def writeResults(gffs, gbks, fas, blasts, outName="termHits.txt"):
    """ Takes an input list for each parameter, and writes each result to the output file """

    with open(outName.name, "w") as out_file:
        if gffs:
            out_file.writelines("==================== GFF3 Term Hits ====================\n\n")
            for gff_hits in gffs:
                out_file.writelines(gff_hits+"\n")
        else:
            pass
        if gbks:
            out_file.writelines("\n==================== GBK Term Hits ====================\n\n")
            for gbk_hits in gbks:
                out_file.writelines(gbk_hits)
        else:
            pass
        if fas:
            out_file.writelines("==================== FASTA Term Hits ====================\n\n")
            for fa_hits in fas:
                out_file.writelines(fa_hits+"\n")
        else:
            pass
        if blasts:
            out_file.writelines("==================== BLAST Term Hits ====================\n\n")
            for blast_hits in blasts:
                out_file.writelines(blast_hits+"\n")
        else:
            pass
            


if __name__ == "__main__":
    print(os.getcwd())
    parser = argparse.ArgumentParser(description="Uses a selection of terms to query an input file for matching cases")
    parser.add_argument("--dbaseTerms",nargs="*",help="dbase terms to search") # will be a select option, based on KEY within the JSON dbase
    parser.add_argument("--custom_txt",nargs="*",help="custom user input terms")
    parser.add_argument("--custom_file",type=argparse.FileType("r"),help="custom new line separated search term file")
    parser.add_argument("--gff3_files",type=argparse.FileType("r"),nargs="*",help="GFF3 File(s)")
    parser.add_argument("--gbk_files",type=argparse.FileType("r"),nargs="*",help="GBK File(s)")
    parser.add_argument("--fa_files",type=argparse.FileType("r"),nargs="*",help="FASTA File(s)")
    parser.add_argument("--blast_files",type=argparse.FileType("r"),nargs="*",help="BLAST.xml File(s)")
    parser.add_argument("--output",type=argparse.FileType("w"),default="termHits.txt")
    args = parser.parse_args()

    ############ STEP I
    ##### Determine user's terms to query
    dbase_terms = dbaseTerms(terms=args.dbaseTerms,galaxy=True)
    user_terms = userTerms(file=args.custom_file,text=args.custom_txt)
    glued_terms = glueTerms(dbase_terms=dbase_terms, user_terms=user_terms)

    ############ STEP II
    ##### Create list with matches
    files = glueFiles(gff=args.gff3_files,gbk=args.gbk_files, fa=args.fa_files, blast=args.blast_files)
    gffs = readGFF3(files=files[0],search_list=glued_terms)
    gbks = readGBK(files=files[1],search_list=glued_terms)
    fas = readFASTA(files=files[2],search_list=glued_terms)
    blasts = readBLAST(files=files[3],search_list=glued_terms)

    ############ STEP III
    ##### Output results to a text file
    writeResults(gffs,gbks,fas,blasts,outName=args.output)


