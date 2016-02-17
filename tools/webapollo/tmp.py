from webapollo import WebApolloInstance, featuresToFeatureSchema
import json
from Bio.SeqFeature import FeatureLocation
from BCBio import GFF

wa = WebApolloInstance('https://cpt.tamu.edu/apollo',
                       'hxr@tamu.edu', 'P@55w0rd')

# fl1 = FeatureLocation(9544, 10484, strand=1)
# fl2 = FeatureLocation(9000, 10484, strand=1)


# for record in wa.bio.ParseRecord('esr.phi29.1'):
    # for feature in record.features:
        # if feature.location.start == 9544:
            # feature.location = fl1
        # elif feature.location.start == 9000:
            # feature.location = fl2

for rec in GFF.parse(open('tmp.gff3', 'r')):
    import pprint; pprint.pprint(featuresToFeatureSchema(rec.features))

# 'annotations', 'dbxrefs', 'description', 'features', 'format', 'id',
# 'letter_annotations', 'lower', 'name', 'reverse_complement', 'seq',
# 'upper'


# wa.annotations.setSequence("Merlin", "0")
# print json.dumps(
    # wa.annotations.getFeatures(),
    # indent=2
# )

# print wa.annotations.getGff3(["303f5cab-4604-483a-a5a8-28fdfd17bfe3"])
# wa.organisms._requestArgs = {'proxies': {'http': 'http://localhost:1000'}}

# import pprint; pprint.pprint(
# wa.organisms.addOrganism('test', '/data/merlin-gff/data')
# )
# print wa.organisms.findAllOrganisms()
# wa.users.createUser('me@me.com', 'firstName', 'lastName', 'password', role='admin')

# oth = wa.users.loadUserById(27)
# wa.users.deleteUser(oth)

# print wa.organisms.deleteOrganism(44)
# print wa.organisms.getSequencesForOrganism('Bob')
# import pprint; pprint.pprint(wa.metrics.getServerMetrics())

# print wa.users.getOrganismPermissionsForUser(15)

# for u in wa.users.loadUsers():
    # # if not u.isAdmin():
        # # wa.users.deleteUser(u)
    # print u, u.isAdmin()
