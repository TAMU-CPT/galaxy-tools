# Remove duplicate terms. As well as add any that is needed.

import explodeJSON as ej
from explodeJSON import save_dict_to_json

lysis_json = "data/lysis-family-expanded.json" # insert json of choice
db = ej.explodeJSON(lysis_json)
db = db.readJSON()

### Add values to dbase:
def add_value_to_term(index_val, db=db, add_value=[]):
    """ index value, put in value """
    for val in add_value:
        db[index_val].append(val)
    
    return db

### Remove values from dbase:
def remove_value_from_term(index_val, db=db, remove_value=[]):
    """ remove values from list """
    for val in remove_value:
        db[index_val].remove(val)
    
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

    db = add_value_to_term("holins", add_value=["anti-holin","anti-holins","co-holin","co-holins"])
    db = remove_value_from_term("holins", remove_value=["lysis (protein)"])
    db = remove_value_from_term("ECDs", remove_value=["viral life cycle","CHAP"])
    save_dict_to_json(obj=db,filename="lysis-family-expanded_culled.json")


