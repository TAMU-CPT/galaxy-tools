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
        db_path = "/galaxy/tools/cpt2/galaxy-tools/tools/proximity/data/lysis-family-v1.0.2.json"
    else:
        #db_path = "/home/adminuser/research/Galaxy-Tools/tools/proximity/data/lysis-family-v1.0.2.json"
        db_path = "data/lysis-family-v1.0.2.json"
    db = ej.explodeJSON(db_path)
    db = db.readJSON()
    dbase_terms = []
    if terms:
        for term in terms:
            index_list = term.split(",")
            for t in index_list:
                if t != "None":
                    dbase_terms.extend(db[t])
                else:
                    dbase_terms = []
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
        if re.search(("__cn__"),str(text[0])):
            #s = text[0].split("__cn__")
            #print(s)
            #print(text[0])
            s = text[0]
            #print(type(s))
            split = s.split("__cn__")
            #print(split)
            user_terms.extend(split)
        else:
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
    """ glue files into one list...I think this is a decent way to go about this...#CHECK LATER#... """
    files = []
    gffs = []
    gbks = []
    blasts = []
    if gff:
        for gff_file in gff:
            gffs.extend(gff_file)
    else:
        pass
    if gbk:
        for gbk_file in gbk:
            gbks.extend(gbk_file)
        #print(gbks)
    else:
        pass
    fas = []
    if fa:
        for fa_file in fa:
            fas.extend(fa_file)
    else:
        pass
    if blast:
        for blast_file in blast:
            blasts.extend(blast_file)
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
        gff3_matches.sort()
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
                    try:
                        if searchInput(str(feature.qualifiers["product"]),search_list=search_list) or searchInput(str(feature.qualifiers["note"]),search_list=search_list) or searchInput(str(feature.qualifiers["dbxref"]),search_list=search_list):
                            gbk_matches.extend([str(feature)])
                        else:
                            continue
                    except KeyError:
                        continue
                gbk_matches = list(set(gbk_matches))
            else:
                print("Parsing - "+file.name)
                record = SeqIO.read(file.name, "genbank")
                for feature in record.features:
                    try:
                        if searchInput(str(feature.qualifiers["product"]),search_list=search_list) or searchInput(str(feature.qualifiers["note"]),search_list=search_list) or searchInput(str(feature.qualifiers["dbxref"]),search_list=search_list):
                            gbk_matches.extend([str(feature)])
                        else:
                            continue
                    except KeyError:
                        continue
                gbk_matches = list(set(gbk_matches))
        gbk_matches.sort()
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
        fa_matches.sort()
        return fa_matches
    else:
        pass

def readBLAST(files,search_list):
    if files:
        for idx, file in enumerate(files):
            if idx == 0:
                print("Parsing - "+file.name)
                blast_records = NCBIXML.parse(open(file.name))
                blast_matches = []
                for blast_record in blast_records:
                    for desc in blast_record.descriptions:
                        pretty = prettifyXML(str(desc))
                        for each_ret in pretty:
                            blast_matches.extend(searchInput(each_ret,search_list=search_list,blast=True,q_id=blast_record.query))
                blast_matches = list(set(blast_matches))
            else:
                print("Parsing - "+file.name)
                blast_records = NCBIXML.parse(open(file.name))
                for blast_record in blast_records:
                    for desc in blast_record.descriptions:
                        pretty = prettifyXML(str(desc))
                        blast_matches.extend(searchInput(each_ret,search_list=search_list,blast=True,q_id=blast_record.query))
                blast_matches = list(set(blast_matches))
            blast_matches.sort()
            return blast_matches
    else:
        pass


######## SEARCH FILE FUNCTIONS
def searchInput(input, search_list,blast=False,q_id=None):
    """ Takes an input search string, and returns uniques of passing """
    output = []
    for search_term in search_list:
        if blast:
            if re.search(re.escape(search_term), input):
                add_query = "QueryID: "+str(q_id)+"\nSearchQuery: "+search_term+"\nMatch: "+input+"\n"
                output.extend([add_query])
            else:
                continue
        #print(search_term)
        #st = r"\b"+search_term+r"\b"
        else:
            if re.search(re.escape(search_term), input):
                #print(search_term+" -> was found")
                output.extend([input])
            else:
                continue
    return list(set(output))

######## prettify-XML function
def prettifyXML(input):
    """ prettifies a string input from a BLAST-xml """
    s = input
    split = s.split(">")

    return split

########## Output File Writer
def writeResults(gffs, gbks, fas, blasts, outName="termHits.txt"):
    """ Takes an input list for each parameter, and writes each result to the output file """

    with open(outName.name, "w+") as out_file:
        if gffs:

            out_file.writelines("\n==================== GFF3 Term Hits ====================\n\n")
            for gff_hits in gffs:
                out_file.writelines(gff_hits+"\n")
        else:
            gffs = []
        if gbks:
            out_file.writelines("\n==================== GBK Term Hits ====================\n\n")
            for gbk_hits in gbks:
                out_file.writelines(gbk_hits+"\n")
        else:
            gbks = []
        if fas:

            out_file.writelines("\n==================== FASTA Term Hits ====================\n\n")
            for fa_hits in fas:
                out_file.writelines(fa_hits+"\n")
        else:
            fas = []
        if blasts:

            out_file.writelines("\n==================== BLAST Term Hits ====================\n\n")
            for blast_hits in blasts:
                out_file.writelines(blast_hits+"\n")
        else:
            blasts = []
        if len(gffs) or len(gbks) or len(fas) or len(blasts):
                print("Terms Found")
        else:
            out_file.writelines("No query matches, try again with new terms!")
            print("No query matches, try again with new terms!")

def write_gff3(gffs,outName="proxHits.gff3"):
    """ writes output to gff3 file for prox2lysis pipeline """

    with open(outName.name, "w+") as out_file:
        out_file.writelines("##gff-version 3\n")
        if gffs:
            for gff_hits in gffs:
                out_file.writelines(gff_hits+"\n")
        else:
            raise Exception("No terms were found from query set")


if __name__ == "__main__":
    print(os.getcwd())
    parser = argparse.ArgumentParser(description="Uses a selection of terms to query an input file for matching cases")
    parser.add_argument("--dbaseTerms",nargs="*",help="dbase terms to search") # will be a select option, based on KEY within the JSON dbase
    parser.add_argument("--custom_txt",nargs="*",help="custom user input terms, if using Galaxy, terms will be __cn__ sep, otherwise by space")
    parser.add_argument("--custom_file",type=argparse.FileType("r"),help="custom new line separated search term file")
    parser.add_argument("--gff3_files",type=argparse.FileType("r"),nargs="*",action="append",help="GFF3 File(s), if multiple files, use another flag")
    parser.add_argument("--gbk_files",type=argparse.FileType("r"),nargs="*",action="append",help="GBK File(s), if multiple files, use another flag")
    parser.add_argument("--fa_files",type=argparse.FileType("r"),nargs="*",action="append",help="FASTA File(s), if multiple files, use another flag")
    parser.add_argument("--blast_files",type=argparse.FileType("r"),nargs="*",action="append",help="BLAST.xml File(s), if multiple files, use another flag")
    parser.add_argument("--output",type=argparse.FileType("w+"),default="termHits.txt")
    parser.add_argument("--prox",action="store_true",help="Use when running the prox2lysis pipeline")
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
    ##### Output results to a text file or gff3
    if args.prox:
        write_gff3(gffs,outName=args.output)
    else:
        writeResults(gffs,gbks,fas,blasts,outName=args.output)


