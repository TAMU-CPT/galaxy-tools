#!/usr/bin/env python
import fileinput

def parse_tsv(trna):
	cols_tsv = trna.split('\t')
	attributes = 'ID="%s";Anticodon="%s";Codon="%s"' % ('tRNA-%s' % cols_tsv[1], cols_tsv[5].lower(), cols_tsv[4])
	cols_gff = [
		cols_tsv[0],
		'aragorn',
		'tRNA',
		cols_tsv[2],
		cols_tsv[3],
		cols_tsv[8],
		'.',
		'.',
		'ID="tRNA-%s";Anticodon="%s";Codon="%s"' % (cols_tsv[1], cols_tsv[5].lower(), cols_tsv[4])
	]
	return cols_gff

if __name__ == '__main__':
	# print file heading
	print '##gff-version-3'
	# process each trna in tsv file
	for trna in fileinput.input():
		print '\t'.join(parse_tsv(trna))