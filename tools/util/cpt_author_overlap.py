#!/usr/bin/env python
import argparse
import logging
logging.basicConfig(level=logging.INFO)


class Paper(object):

    def __init__(self, id, authors):
        self.id = id
        self.authors = self.process(authors)

    def process(self, author_list):
        full_list = []
        if ' and ' in author_list:
            items = author_list.split(' and ')
            full_list.append(items[1])
            full_list.extend(items[0].split(', '))
        else:
            full_list = author_list.split(', ')  # just in case...
        return full_list

    def overlap(self, other_paper):
        count = 0
        for author in self.authors:
            for other_author in other_paper.authors:
                if author == other_author:
                    count += 1
        return count

    def close_enough(self, other_paper, threshold=3):
        return self.overlap(other_paper) >= threshold

    def __str__(self):
        return '%s by %s' % (self.id, ', '.join(self.authors))


def overlap(table):
    groups = [
    ]
    for line in table.readlines():
        if not line.startswith('#'):
            linedata = line.strip().split('\t')
            paper = Paper(linedata[0], linedata[1])
            added = False
            for group in groups:
                if True in [paper.close_enough(x) for x in group]:
                    group.append(paper)
                    added = True
                    continue
            if not added:
                groups.append([paper])

    ids = []
    for group in groups:
        ids.append([x.id for x in group])
    return ids

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Determine overlaps in CSV in colum data')
    parser.add_argument('table', type=file, help='2 column table')
    parser.add_argument('--version', action='version', version='0.1')
    args = parser.parse_args()

    data = overlap(**vars(args))
    for group in data:
        print '\t'.join(group)
