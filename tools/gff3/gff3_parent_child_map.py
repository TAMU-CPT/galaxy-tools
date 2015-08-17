import pprint
import sys
from BCBio.GFF import GFFExaminer

with open(sys.argv[1], 'r') as handle:
    examiner = GFFExaminer()
    pprint.pprint(examiner.parent_child_map(handle))
