# Proximity (to Lysis) Scripts

## `generateLysisFamily.py`
* Has numerous list with lysis family relationships. Output is a dictionary --> json that stores these relationships for future use.

## `synonymParse.py`
* module to incooperate quick parsing of new-line separated files
```python
import synonymParse as sp
file = 'path/to/file.txt'

myObject = sp.Synonym(filename = file, delims = '\t)
myObjectList = sp.myObject.parse_it() # creates a list of items from the input file

# Can iterate across the list, and will be able to do various manipulations with it
```