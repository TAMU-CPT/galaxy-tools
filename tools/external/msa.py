#!/usr/bin/env python
import sys
import argparse


class MSA:
    """
    multiple item alignment
    """

    def __init__(self, bidi=True, gap=0, match=5, mismatch=-1):
        self.sequences = []
        self.bidi = bidi
        self.relationships = {}
        self.gap_penalty = gap
        self.match_score = match
        self.mismatch_score = mismatch

        self.merger = []
        self.number_of_aligned_lists = 0

    def add_relationship(self, a, b):
        if a not in self.relationships:
            self.relationships[a] = []
        self.relationships[a].append(b)

        if self.bidi:
            if b not in self.relationships:
                self.relationships[b] = []
            self.relationships[b].append(a)

    def Sij(self, merger_row, query):
        for elem in merger_row:
            if (query in self.relationships and elem in self.relationships[query]) or (
                elem in self.relationships and query in self.relationships[elem]
            ):
                return self.match_score
        return self.mismatch_score

    def align_list(self, data):
        # If we haven't aligned any lists, we do something special
        if self.number_of_aligned_lists == 0:
            # Pretend we've aligned ONE list already
            self.number_of_aligned_lists = 1
            self.merger = [[x] for x in data]
        else:
            self.find_best_path(data)
            self.number_of_aligned_lists += 1

    def find_best_path(self, data):
        max_i = len(self.merger)
        max_j = len(data)

        score_mat = {}
        point_mat = {}

        point_mat[(0, 0)] = "Z"
        score_mat[(0, 0)] = 0

        for i in range(max_i):
            point_mat[(i + 1, 0)] = "U"
            score_mat[(i + 1, 0)] = self.gap_penalty

        for i in range(max_j):
            point_mat[(0, i + 1)] = "L"
            score_mat[(0, i + 1)] = self.gap_penalty

        # Score
        for i in range(max_i):
            ci = self.merger[i]
            for j in range(max_j):
                cj = data[j]
                # scoring
                diag_score = score_mat[(i, j)] + self.Sij(ci, cj)
                up_score = score_mat[(i, j + 1)] + self.gap_penalty
                left_score = score_mat[(i + 1, j)] + self.gap_penalty

                pos = (i + 1, j + 1)
                if diag_score >= up_score:
                    if diag_score >= left_score:
                        score_mat[pos] = diag_score
                        point_mat[pos] = "D"
                    else:
                        score_mat[pos] = left_score
                        point_mat[pos] = "L"
                else:
                    if up_score >= left_score:
                        score_mat[pos] = up_score
                        point_mat[pos] = "U"
                    else:
                        score_mat[pos] = left_score
                        point_mat[pos] = "L"

        # self.print2dArray(score_mat, data)
        # self.print2dArray(point_mat, data)

        new_row_set = []
        i = max_i + 0
        j = max_j + 0
        while True:
            if i == 0 and j == 0:
                break

            d = point_mat[(i, j)]
            new_row = None
            if d == "D":
                new_row = self.merger[i - 1] + [data[j - 1]]
                i -= 1
                j -= 1
            elif d == "L":
                new_row = ["-" for _ in range(self.number_of_aligned_lists)]
                new_row.append(data[j - 1])
                j -= 1
            elif d == "U":
                new_row = self.merger[i - 1] + ["-"]
                i -= 1

            new_row_set.append(new_row)

        self.merger = new_row_set[::-1]

    def print2dArray(self, d, xlab=None):
        k = d.keys()
        i = max([q[0] for q in k])
        j = max([q[1] for q in k])
        if xlab:
            sys.stderr.write("  ".join(["x"] + [str(z) for z in xlab]))
        else:
            sys.stderr.write("  ".join(["x"] + [str(z) for z in range(j)]))
        sys.stderr.write("\n")
        for n in range(i):
            sys.stderr.write("%s] " % n)
            for m in range(j):
                sys.stderr.write(str(d.get((n, m), "?")) + " ")
            sys.stderr.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rels", type=argparse.FileType("r"), help="Relationships")
    parser.add_argument("lists", type=argparse.FileType("r"), nargs="+", help="Lists")
    ARGS = parser.parse_args()

    m = MSA()
    for line in ARGS.rels:
        d = line.strip().split()
        m.add_relationship(*d)

    for l in ARGS.lists:
        m.align_list([x.strip() for x in l])

    # Done aligning
    for row in m.merger:
        sys.stdout.write("\t".join(row))
        sys.stdout.write("\n")

    # Test Case
    # m = MSA()
    # m.add_relationship('1', 'a')
    # m.add_relationship('2', 'b')
    # m.add_relationship('3', 'd')

    # m.add_relationship('1', 'x')
    # m.add_relationship('2', 'y')
    # m.add_relationship('2', 'z')

    # m.align_list(list('abcdef'))
    # m.align_list(list('12345'))
    # for row in m.merger:
    # print '\t'.join(row)

    # print
    # m.align_list(list('uvwxyz'))
