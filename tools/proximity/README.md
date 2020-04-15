# Proximity (to Lysis) Scripts
## `generateLysisFamily.py`
* Numerous list with lysis family relationships. Output is a dictionary --to-a--> json that stores these relationships for future use.

## `explodeJSON.py`
* Quick and easy class that will permit (easy) json manipulation.

## `fileProx.py`
* GALAXY tool that allows users to query an input file for a set of search terms; either custom, dbase, or a combination of the two.

## `goQuery.py`
* Takes terms from the lysis-family.json and queries quickgo

## `goSynonymAdd.py`
* Merges the lysis-family.json object with results from quickgo

# `/data`
* Check README of this directory to see details on the `.txt` and `.json` files.

# `/test-data`
* test.json for testing various scripts. Used for `goQuery.py`
* sample.gff3
* sample.fa
* sample.xml
* search.txt