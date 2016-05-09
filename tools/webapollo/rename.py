from webapollo import WebApolloInstance, featuresToFeatureSchema
import json
from Bio.SeqFeature import FeatureLocation

wa = WebApolloInstance('https://cpt.tamu.edu/apollo',
                       'hxr@tamu.edu', 'P@55w0rd')

wa.annotations.setSequence("DUC3_BJones", "DUC3_BJones")

with open('duc3-bjones.gff3', 'r') as handle:
    for idx, row in enumerate(handle):
        (a, name, b) = row.strip().split('\t')
        # wa.annotations.setName(a, name)
        # wa.annotations.setName(b, name)
        try:
            wa.annotations.setName(a, name)
        except:
            print "Failed on " + a

        try:
            wa.annotations.setName(b, name)
        except:
            print "Failed on " + b
        print idx


# print wa.annotations.setName("3c0c2904-7643-4b48-b57f-fe1d24771340", "TEST")
# print wa.annotations.setName("ed32f9b3-6f1e-4bb3-b424-2cbdc4e7e4a0", "TEST")


# # print json.dumps(
    # # wa.annotations.getFeatures(),
    # # indent=2
# # )
