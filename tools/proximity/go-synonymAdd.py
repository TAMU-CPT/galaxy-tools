########## parses a result from the goQuery.py and appends results and syns to data

import json
import pandas as pd 
import explodeJSON as ej
import re


dbase = "data/lysis-family.json" # generated from generateLysisFamily.py
go_file = "go-synonym-results-JRR.txt" # cleared up results from JRR

### Using ej to parse json
db = ej.explodeJSON(dbase)
db = db.readJSON()

### Using pandas to peek into txt file
cols = ["query_term", "GO:id", "GO:name", "description", "restriction", "synonyms"]
go = pd.read_csv(go_file, names=cols, sep="\t")
go_syns = go[["query_term","GO:name","synonyms"]] # subset of what I want from the datafile

def parse_pipes(piped_data): # PASS
    """ 
    Extracts synonyms, and separates them, based on the pipe delimiter used from goQuery.py
    """
    pipe = piped_data.split("|")
    return pipe

def remove_fat(term): # PASS
    """
    removes extra fluff of terms, such as regulation, etc... 
    Do this within a loop to suck all the fat out of terms/synonyms
    """
    s = term
    if re.search(("regulation"), s):
        s = s.split("regulation of")[1]
        s = s.strip()
    else:
        pass
    if re.search(("activity"), s):
        s = s.split(" activity")[0]
        s = s.strip()
    else:
        pass
    if re.search(("inhibition of"), s):
        s = s.split("inhibition of")[1]
        s = s.strip()
    else:
        pass
    if re.search(("None"), s):
        s = ""
    else:
        pass
    if re.search(("NoneAtAll"), s):
        s = ""
    else:
        pass
    return s

def readFrameVals(frame=go_syns): # PASS
    """
    reads the go_syns frame, and grabs the values we want from it
    """
    link = frame["query_term"].tolist()
    go = frame["GO:name"].tolist()
    syn = frame["synonyms"].tolist()
    return link, go, syn, zip(link,go,syn)

def form_dict(link, go, syn): # PASS
    """
    format a dictionary that will stack all matching query terms (eventually mapped back to json); remove duplicate / redundant names and synonyms
    """
    sync = {}
    for i, name in enumerate(link):
        try:
            if sync[name] != []:
                sync[name].append(remove_fat(go[i]))
                pipes = parse_pipes(syn[i])
                for p in pipes:
                    p = remove_fat(p)
                    sync[name].append(p)
        except KeyError:
            sync[name] = []
            sync[name].append(remove_fat(go[i]))
            pipes = parse_pipes(syn[i])
            for p in pipes:
                p = remove_fat(p)
                sync[name].append(p)
    # Removing Duplicates
    unq = {}
    for k, v in sync.items():
        unq[k] = list(set(v))
        # Removing "None" and "NoneAtAll", this should be the last cleaning needed
        is_not = lambda x: x is not "None" or "NoneAtAll"
        unq[k] = list(filter(is_not, unq[k]))

    return unq

def match_and_store(input, db=db): # FAIL
    """ return JSON object that has go-synonyms added ; input will be return of form_dict function"""
    """
    for mapper, synonyms in input.items():
        for matching_mapper_list in dbase.values():
            for matching_mapper in matching_mapper_list:
                if matching_mapper == mapper:
                    matching_mapper_list += synonyms
                else:
                    continue
    """
    new_db = {}
    for mapper, synonyms in input.items():
        #print(mapper)
        #print("//////////////////")
        #print(synonyms)
        for family_term, list_of_syns in db.items():
            new_db[family_term] = [] # possibly format db
            #print(family_term)
            #print("++++++++++++++++")
            #print(list_of_syns)
            for map_synonym in list_of_syns:
                #print(map_synonym)
                if mapper == map_synonym:
                    print(mapper+" should equal "+map_synonym)
                    #print("should equal "+map_synonym)
                    #print(map_synonym)
                    print(synonyms)
                    list_of_syns.append(synonyms)
                    #list_of_syns += synonyms
                else:
                    continue
    #print(dbase)
    """
    for vals in db.values():
        print(vals)
        for v in vals:
            print(v)
    """
    return db

if __name__ == "__main__":
    #term = "up-regulation of mucopeptide N-acetylmuramoylhydrolase activity"
    #pipes = "up-regulation of mucopeptide N-acetylmuramoylhydrolase activity|positive regulation of 1,4-N-acetylmuramidase activity|up regulation of muramidase activity|up-regulation of mucopeptide glucohydrolase activity|up-regulation of muramidase activity|up-regulation of peptidoglycan N-acetylmuramoylhydrolase activity|up regulation of peptidoglycan N-acetylmuramoylhydrolase activity|upregulation of mucopeptide N-acetylmuramoylhydrolase activity|upregulation of mucopeptide glucohydrolase activity|upregulation of peptidoglycan N-acetylmuramoylhydrolase activity|positive regulation of mucopeptide glucohydrolase activity|positive regulation of muramidase activity|up regulation of mucopeptide N-acetylmuramoylhydrolase activity|up regulation of mucopeptide glucohydrolase activity|up-regulation of lysozyme activity|upregulation of lysozyme activity|positive regulation of N,O-diacetylmuramidase activity|positive regulation of mucopeptide N-acetylmuramoylhydrolase activity|up-regulation of 1,4-N-acetylmuramidase activity|up-regulation of N,O-diacetylmuramidase activity|upregulation of N,O-diacetylmuramidase activity|positive regulation of peptidoglycan N-acetylmuramoylhydrolase activity|upregulation of 1,4-N-acetylmuramidase activity|upregulation of muramidase activity|up regulation of 1,4-N-acetylmuramidase activity|up regulation of N,O-diacetylmuramidase activity|up regulation of lysozyme activity"
    #f = remove_fat(term)
    #print(f)
    #vals = parse_pipes(pipes)
    #print(vals)
    #for v in vals:
    #    remove_fat(v)
    link, go, syn, zipped = readFrameVals()
    #print(syn)
    u = form_dict(link=link,go=go,syn=syn)
    
    #print(u)

    #print(type(db))

    #print(db)
    
    #for k, v in db.items():
    #    print(k)
    #    print(len(v))

    complete = match_and_store(input=u,db=db)

    #for k,v in complete.items():
    #    print(k)
    #    print(v)

    # SAVE AS JSON
    #filename = "lysis-family-expanded.json"
    #with open("data/" + filename, "w") as j:
    #    json.dump(complete, j, indent="\t")
