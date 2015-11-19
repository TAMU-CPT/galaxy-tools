from webapollo import WebApolloInstance
# from webapollo import UserObj
import random

wa = WebApolloInstance('http://localhost:8080/apollo',
                       'hxr@tamu.edu', 'pass')
r = random.randint(1, 1000)

me = wa.users.loadUserById(15)
# print me
wa.users.updateUser(me, 'hxr@tamu.edu', 'Helena', 'Rasche', 'pass')
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
