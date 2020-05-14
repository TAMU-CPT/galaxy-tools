# Proximity (to Lysis) Scripts
_aside: Some of the scripts that are prepatory scripts, have been moved, and thus their paths for some of the calls they make will raise errors. When/If this is ever ran from the ground up again, that will need to be remedied first_

_aside#2: planemo test will need to be regenerated because of the XML alterations, and need to test the prox function. I'll still likely just stick with one test case_

### `searchFile.py`
* Script that searches a multitude of filetypes, and multiple at a time, and queries the file(s) for terms selected by the user.
* `--prox` flag enables the proximity2lysis GFF3 output

### `editDB.py`
* Script that finds redudancy, and other issues, with the `...expanded.json` output and resolves said issues.
* Functions to modify the database, and output a new version.

## `/prep-data`
### `generateLysisFamily.py`
* Numerous list with lysis family relationships. Output is a dictionary --to-a--> json that stores these relationships for future use.

### `explodeJSON.py`
* Quick and easy class that will permit (easy) json manipulation.

### `goQuery.py`
* Takes terms from the lysis-family.json and queries quickgo

### `goSynonymAdd.py`
* Merges the lysis-family.json object with results from quickgo

## `/data`
* Check README of this directory to see details on the `.txt` and `.json` files.
* Order of scripts to generate `lysis-family-expanded.json` (this is the dbase used in query scripts/galaxy tools)
    * `generateLysisFamily.py`
    * `goQuery.py` + curate results --- _then_ ---> `goSynonymAdd.py`
    * I don't feel like the list was appropriately parsed down, especially during the merger. Thus run the output json from the above through `editDB.py`

## `/test-data`
* `/other`
    * test.json for testing various scripts. Used for `goQuery.py`
    * sample.gff3
    * sample.fa
    * sample.xml
* `/general_use` :: for running planemo test as well as testing the overall function of the script
* `/prox` :: for testing the proximity to lysis portion of this tool

## Order for reproducing database
1. generateLysisFamily.py
2. goQuery.py
3. goSynonymAdd.py
4. (optional) editDB.py

# Release Notes:
CPT: OTHER- Search File
In light of various lysis projects, searching for lysis-related terms within files can be useful. This tool allows for the querying of GFF3, FASTA, Genbank, and/or BLAST-xml files with a user defined set of terms. The terms can either be from the curated lysis family term database, a custom set of terms, or both. Each query term is passed to each file, to see if it is within the contents (similar to CTRL+F), but it queries specific places, specific to the file type. For example, genbank’s product and notes qualifier is queried. Matches of each filetype are returned in an output file that allows users to see the contents of what was matched.
