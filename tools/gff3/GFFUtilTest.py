#!/usr/bin/env python
import sys
import gffutils

db = gffutils.create_db('test-data/cds.gff3', dbfn='temp.db', force=True, keep_order=True,merge_strategy='merge', sort_attribute_values=True)
db2 = gffutils.FeatureDB('temp.db', keep_order=True)
#print(dir(db2))
#print((db2.directives))
#gene = db2['cds_orf00001']
#print(dir(gene))
for y in db2.directives:
  print('##' + y)
for x in db2.all_features():
  print(x)
