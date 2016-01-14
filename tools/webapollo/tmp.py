from webapollo import WebApolloInstance
import json

wa = WebApolloInstance('http://localhost:8080/apollo-2.0.2-SNAPSHOT',
                       'hxr@tamu.edu', 'p')


wa.annotations.setSequence("Merlin", "0")
print json.dumps(
    wa.annotations.getFeatures(),
    indent=2
)

print wa.annotations.getGff3(["303f5cab-4604-483a-a5a8-28fdfd17bfe3"])
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
