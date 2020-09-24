# Remove duplicate terms. As well as add any that is needed.

import explodeJSON as ej
from explodeJSON import save_dict_to_json


### create new key
def add_new_key(db, add_key=[]):
    """ Set of keys to add to the database """
    for new_key in add_key:
        db[new_key] = []

    return db


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

    lysis_json = "data/lysis-family-v1.0.2.json" # insert json of choice
    db = ej.explodeJSON(lysis_json)
    db = db.readJSON()
    #revise_db = add_new_key(db=db,add_key=["spanins"])
    #files = ["data/term_additions/200505_holin_domains.txt","data/term_additions/200505_Spanin_Domains.txt"]
    terms = ["DUF2570","PF10828","IPR022538","DUF2514","PF10721","IPR019659","DUF2681","PF10883","IPR020274"]
    #revise_db = add_from_file(files[0],"holin_domains",revise_db)
    #revise_db = add_from_file(files[1],"spanin_domains",revise_db)
    revise_db = add_value_to_term("spanin_domains",db,add_value=terms)
    save_dict_to_json(obj=revise_db,filename="data/lysis-family-v1.0.3.json")


