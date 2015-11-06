import requests
import json

class WebApolloInstance(object):

    def __init__(self, url, username, password):
        self.apollo_url = url
        self.username = username
        self.password = password

        self.annotations = AnnotationsClient(self)
        self.groups = GroupsClient(self)
        self.io = IOClient(self)
        self.organisms = OrganismsClient(self)
        self.users = UsersClient(self)


class GroupObj(object):
    def __init__(self, **kwargs):
        self.name = kwargs['name']


class UserObj(object):
    ROLE_USER = 'USER'
    ROLE_ADMIN = 'ADMIN'

    def __init__(self, **kwargs):
        for attr in kwargs.keys():
            setattr(self, attr, kwargs[attr])

        if 'groups' in kwargs:
            groups = []
            for groupData in kwargs['groups']:
                groups.append(GroupObj(**groupData))
            self.groups = groups

        self.__props = kwargs.keys()


    def isAdmin(self):
        if hasattr(self, 'role'):
            return self.role == self.ROLE_ADMIN
        return False

    def refresh(self, wa):
        # This method requires some sleeping usually.
        newU = wa.users.loadUser(self).toDict()
        for prop in newU:
            setattr(self, prop, newU[prop])

    def toDict(self):
        data = {}
        for prop in self.__props:
            data[prop] = getattr(self, prop)
        return data

    def __str__(self):
        return '<User %s: %s %s <%s>>' % (self.userId, self.firstName,
                                          self.lastName, self.username)


class Client(object):

    def __init__(self, webapolloinstance, **requestArgs):
        self.wa = webapolloinstance

        self.verify = requestArgs.get('verify', True)
        self.requestArgs = requestArgs

        if 'verify' in self.requestArgs:
            del self.requestArgs['verify']

    def request(self, clientMethod, data, post_params={}):
        url = self.wa.apollo_url + self.CLIENT_BASE + clientMethod

        headers = {
            # 'Content-Type': 'application/json'
            'Content-Type': 'application/x-www-form-urlencoded',
        }

        data.update({
            'username': self.wa.username,
            'password': self.wa.password,
        })

        formatted_data = "data=%s" % json.dumps(data)
        r = requests.post(url, data=formatted_data, headers=headers,
                          verify=self.verify, params=post_params)
        if r.status_code == 200:
            return r.json()
        # @see self.body for HTTP response body
        raise Exception("Unexpected response from apollo %s: %s" %
                        (r.status_code, r.text))


class AnnotationsClient(Client):
    CLIENT_BASE = '/annotationEditor/'


class GroupsClient(Client):
    CLIENT_BASE = '/group/'


class IOClient(Client):
    CLIENT_BASE = '/ioService/'

    def write(self, exportType='FASTA', seqType='peptide',
              exportFormat='text', sequences=None, organism=None,
              output='text', exportAllSequences=False,
              exportGff3Fasta=False):
        if exportType not in ('FASTA', 'GFF3'):
            raise Exception("exportType must be one of FASTA, GFF3")

        if seqType not in ('peptide', 'cds', 'cdna', 'genomic'):
            raise Exception("seqType must be one of peptide, cds, dna, genomic")

        if exportFormat not in ('gzip', 'text'):
            raise Exception("exportFormat must be one of gzip, text")

        if output not in ('file', 'text'):
            raise Exception("output must be one of file, text")

        data = {
            'type': exportType,
            'seqType': seqType,
            'format': exportFormat,
            'sequences': sequences,
            'organism': organism,
            'output': output,
            'exportAllSequences': exportAllSequences,
            'exportGff3Fasta': exportGff3Fasta,
        }

        return self.request('write', data)

    def download(self, uuid, outputFormat='gzip'):

        if outputFormat.lower() not in ('gzip', 'text'):
            raise Exception("outputFormat must be one of file, text")

        data = {
            'format': outputFormat,
            'uuid': uuid,
        }
        return self.request('write', data)


class OrganismsClient(Client):
    CLIENT_BASE = '/organism/'

    def addOrganism(self, commonName, directory, blatdb=None, species=None, genus=None, public=False):
        data = {
            'commonName': commonName,
            'directory': directory,
            'publicMode': public,
        }

        if blatdb is not None:
            data['blatdb'] = blatdb
        if genus is not None:
            data['genus'] = genus
        if species is not None:
            data['species'] = species

        return self.request('addOrganism', data)

    def findAllOrganisms(self):
        return self.request('findAllOrganisms', {})

    def deleteOrganism(self, organismId):
        return self.request('deleteOrganism', {'id': organismId})

    def deleteOrganismFeatures(self, organismId):
        return self.request('deleteOrganismFeatures', {'id': organismId})

    def getSequencesForOrganism(self, commonName):
        return self.request('getSequencesForOrganism', {'organism': commonName})

    def updateOrganismInfo(self, organismId, commonName, directory, blatdb=None, species=None, genus=None, public=False):
        data = {
            'id': organismId,
            'commonName': commonName,
            'directory': directory,
            'publicMode': public,
        }

        if blatdb is not None:
            data['blatdb'] = blatdb
        if genus is not None:
            data['genus'] = genus
        if species is not None:
            data['species'] = species

        return self.request('updateOrganismInfo', data)


class UsersClient(Client):
    CLIENT_BASE = '/user/'

    def getOrganismPermissionsForUser(self, user):
        data = {
            'userId': user.userId
        }
        data
        # return self.request('getOrganismPermissionsForUser', data)

    def updateOrganismPermission(self, user, organism, administrate=False,
                                 write=False, export=False, read=False):
        data = {
            'userId': user.userId,
            'organism': organism,
            'administrate': administrate,
            'write': write,
            'export': export,
            'read': read,
        }
        data
       # return self.request('updateOrganismPermission', data)

    def loadUser(self, user):
        return self.loadUserById(user.userId)

    def loadUserById(self, userId):
        res = self.request('loadUsers', {'userId': userId})
        if isinstance(res, list):
            # We can only match one, right?
            return UserObj(**res[0])
        else:
            return res

    def loadUsers(self, userId=None):
        data = {}
        if userId is not None:
            data['userId'] = userId

        res = self.request('loadUsers', data)
        if isinstance(res, list):
            return [UserObj(**x) for x in res]
        else:
            return res

    def addUserToGroup(self, group, user):
        data = {'group': group.name, 'userId': user.userId}
        return self.request('addUserToGroup', data)

    def removeUserFromGroup(self, group, user):
        data = {'group': group.name, 'userId': user.userId}
        return self.request('removeUserFromGroup', data)

    def createUser(self, email, firstName, lastName, newPassword, role="user", groups=None):
        data = {
            'firstName': firstName,
            'lastName': lastName,
            'email': email,
            'role': role,
            'groups': [] if groups is None else groups,
            # 'availableGroups': [],
            'newPassword': newPassword,
            # 'organismPermissions': [],
        }
        return self.request('createUser', data)

    def deleteUser(self, user):
        return self.request('deleteUser', {'userId': user.userId})

    def updateUser(self, user, email, firstName, lastName, newPassword):
        data = {
            'userId': user.userId,
            'email': email,
            'firstName': firstName,
            'lastName': lastName,
            'newPassword': newPassword,
        }
        return self.request('updateUser', data)
