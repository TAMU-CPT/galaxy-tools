#!/usr/bin/env python
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='split blast results by organism')
    parser.add_argument('blast', type=argparse.FileType("r"))
    args = parser.parse_args()

    blast = [
        x.split('\t')
        for x in args.blast.readlines()
    ]

    for row in blast:
        if '<>' in row[24]:
            for i in row[24].strip().split('<>'):
                row_copy = row
                row_copy[24] = i
                print '\t'.join(row_copy).strip()
        else:
            print '\t'.join(row).strip()
