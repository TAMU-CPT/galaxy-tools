# Remove duplicate terms. As well as add any that is needed.

import explodeJSON as ej
from explodeJSON import save_dict_to_json


lysis_json = "data/lysis-family-expanded.json" # insert json of choice
db = ej.explodeJSON(lysis_json)
db = db.readJSON()

### Add values to dbase:
def add_value_to_term(index_val, db, add_value=[]):
    """ index value, put in value """
    for val in add_value:
        db[index_val].append(val)
    
    return db

### Remove values from dbase:
def remove_value_from_term(index_val, db, remove_value=[]):
    """ remove values from list """
    for val in remove_value:
        db[index_val].remove(val)
    
    return db

### Terms to add from a file
def add_from_file(input_file,index_val,db,sep="\n"):
    """ input file, new line separated currently, and append files to correct key, return is altered dictionary"""
    terms = open(input_file).read().splitlines()
    db = add_value_to_term(index_val,db,terms)
    return db

if __name__ == "__main__":
    ### Take a peak at numbers
    for k, v in db.items():
        print(k)
        print(len(v))
        print(len(set(v)))
        v = list(set(v))
    ### Replace with the unique set of the internal synonym list
    for k, v in db.items():
        db[k] = list(set(v))
    ### verify lengths (comparing to first iteration)
    for k,v in db.items():
        print("new key --> "+str(k))
        print("new length --> "+str(len(v)))

    db = add_value_to_term("holins", db, add_value=["anti-holin","anti-holins","co-holin","co-holins"])
    db = remove_value_from_term("holins", db, remove_value=["lysis (protein)"])
    db = remove_value_from_term("ECDs", db, remove_value=["viral life cycle","CHAP"])
    db = remove_value_from_term("ECDs_accro", db, remove_value=["DUF"])
    
    # Modify database as per #374 issue
    revise_db = {"endolysins":[],"holins":[],"endolysin_domains":[]}
    for k, v in db.items():
        if k == "endolysins":
            revise_db["endolysins"].extend(v)
        elif k == "holins":
            revise_db["holins"].extend(v)
        else:
            revise_db["endolysin_domains"].extend(v)
    
    ### Add new terms from @LMad
    file = "data/endolysin domains to add.txt"
    revise_db = add_from_file(file,"endolysin_domains",revise_db)
    
    save_dict_to_json(obj=revise_db,filename="data/lysis-family-v1.0.0.json")


