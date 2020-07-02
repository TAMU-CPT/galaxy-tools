#!/usr/bin/env python
import fileinput

print "##gff-version-3"
# process each trna in tsv file
metaLines = 0
for trna in fileinput.input():
    if metaLines < 3:
      metaLines += 1
      continue
    cols_tsv = trna.split("\t")
    cols_gff = [
        cols_tsv[0],
        "tRNAscan-SE",
        "tRNA",
        cols_tsv[2],
        cols_tsv[3],
        cols_tsv[8],
        ".",
        ".",
        'ID="trna.%s";Anticodon="%s";Codon="tRNA-%s"'
        % (cols_tsv[1], cols_tsv[5].lower(), cols_tsv[4]),
    ]
    print "\t".join(cols_gff)
