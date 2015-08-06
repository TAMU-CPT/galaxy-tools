#!/usr/bin/env python
import fileinput
import sys

def parse_tsv(trna):
	cols_gff = ['.' for i in range(9)]
	cols_tsv = trna.split()
	cols_gff[0] = cols_tsv[0]
	cols_gff[1] = 'aragorn'
	cols_gff[2] = 'tRNA'
	cols_gff[3] = cols_tsv[2]
	cols_gff[4] = cols_tsv[3]
	cols_gff[5] = cols_tsv[8]
	attributes = 'ID="tRNA-%s' % cols_tsv[1]
	attributes += '";Anticodon="' + cols_tsv[5].lower() + '";"Codon="' + cols_tsv[4] + '"'
	cols_gff[8] = attributes
	return cols_gff

if __name__ == '__main__':
	# print file heading
	print '##gff-version-3'
	# process each trna in tsv file
	for trna in fileinput.input(sys.argv[1]):
		print '\t'.join(parse_tsv(trna))