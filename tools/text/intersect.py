#!/usr/bin/env python
import argparse
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def joinStripRow(row):
    return "\t".join(row).strip()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("a", type=argparse.FileType("r"))
    parser.add_argument("a_col", type=int, default=1)
    parser.add_argument("b", type=argparse.FileType("r"))
    parser.add_argument("b_col", type=int, default=1)
    parser.add_argument("--invert", action="store_true")
    parser.add_argument("--output_a", type=argparse.FileType("w"), default="out_a.tab")
    parser.add_argument("--output_b", type=argparse.FileType("w"), default="out_b.tab")
    parser.add_argument(
        "--output_ab", type=argparse.FileType("w"), default="out_ab.tab"
    )
    args = parser.parse_args()

    data_a = [x.split("\t") for x in args.a.readlines()]

    data_b = [x.split("\t") for x in args.b.readlines()]

    data_a_indexed = {}
    for row in data_a:
        if row[args.a_col - 1] not in data_a_indexed:
            data_a_indexed[row[args.a_col - 1]] = [row]
        else:
            data_a_indexed[row[args.a_col - 1]].append(row)

    marked_rows_a = []
    marked_rows_b = []

    for num, row in enumerate(data_b):
        row_id = row[args.b_col - 1]
        if row_id in data_a_indexed:

            marked_rows_b.append(num)
            if row_id not in marked_rows_a:
                marked_rows_a.append(row_id)

            if not args.invert:
                for i in data_a_indexed[row_id]:
                    args.output_ab.write(
                        joinStripRow(i) + "\t" + joinStripRow(row) + "\n"
                    )

    if args.invert:
        for i in data_a_indexed:
            if i not in marked_rows_a:
                for j in data_a_indexed[i]:
                    args.output_a.write(joinStripRow(j) + "\n")
        for num, i in enumerate(data_b):
            if num not in marked_rows_b:
                args.output_b.write(joinStripRow(i) + "\n")
