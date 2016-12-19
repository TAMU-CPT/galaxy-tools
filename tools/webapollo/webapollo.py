import requests
import json
import os
import collections
import StringIO
import logging
import argparse
from BCBio import GFF
from Bio import SeqIO
logging.getLogger("requests").setLevel(logging.CRITICAL)
log = logging.getLogger()


def WAAuth(parser):
    parser.add_argument('apollo', help='Complete Apollo URL')
    parser.add_argument('username', help='WA Username')
    parser.add_argument('password', help='WA Password')
    parser.add_argument('--remote_user', default='', help='If set, ignore password, set the header with the name supplied to this argument to the value of email')


def OrgOrGuess(parser):
    parser.add_argument('--org_json', type=argparse.FileType("r"), help='Apollo JSON output, source for common name')
    parser.add_argument('--org_raw', help='Common Name')
    parser.add_argument('--org_id', help='Organism ID')


def CnOrGuess(parser):
    OrgOrGuess(parser)
    parser.add_argument('--seq_fasta', type=argparse.FileType("r"), help='Fasta file, IDs used as sequence sources')
    parser.add_argument('--seq_raw', nargs='*', help='Sequence Names')


def GuessOrg(args, wa):
    if args.org_json:
        orgs = [x.get('commonName', None)
                for x in json.load(args.org_json)]
        orgs = [x for x in orgs if x is not None]
        return orgs
    elif args.org_raw:
        org = args.org_raw.strip()
        if len(org) > 0:
            return [org]
        else:
            raise Exception("Organism Common Name not provided")
    elif args.org_id:
        return [wa.organisms.findOrganismById(args.org_id).get('commonName', None)]
    else:
        raise Exception("Organism Common Name not provided")


def GuessCn(args, wa):
    org = GuessOrg(args, wa)
    seqs = []
    if args.seq_fasta:
        # If we have a fasta, pull all rec ids from that.
        for rec in SeqIO.parse(args.seq_fasta, 'fasta'):
            seqs.append(rec.id)
    elif args.seq_raw:
        # Otherwise raw list.
        seqs = [x.strip() for x in args.seq_raw if len(x.strip()) > 0]

    return org, seqs


def AssertUser(user_list):
    if len(user_list) == 0:
        raise Exception("Unknown user. Please register first")
    else:
        return user_list[0]


def AssertAdmin(user):
    if user.role == 'ADMIN':
        return True
    else:
        raise Exception("User is not an administrator. Permission denied")


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
        self.metrics = MetricsClient(self)
        self.bio = RemoteRecord(self)

    def __str__(self):
        return '<WebApolloInstance at %s>' % self.apollo_url


class GroupObj(object):

    def __init__(self, **kwargs):
        self.name = kwargs['name']

        if 'id' in kwargs:
            self.groupId = kwargs['id']


class UserObj(object):
    ROLE_USER = 'USER'
    ROLE_ADMIN = 'ADMIN'

    def __init__(self, **kwargs):
        # Generally expect 'userId', 'firstName', 'lastName', 'username' (email)
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
        self._wa = webapolloinstance

        self.__verify = requestArgs.get('verify', True)
        self._requestArgs = requestArgs

        if 'verify' in self._requestArgs:
            del self._requestArgs['verify']

    def request(self, clientMethod, data, post_params={}, isJson=True):
        url = self._wa.apollo_url + self.CLIENT_BASE + clientMethod

        headers = {
            'Content-Type': 'application/json'
        }

        data.update({
            'username': self._wa.username,
            'password': self._wa.password,
        })

        r = requests.pist(url, data=json.dumps(data), headers=headers,
                          verify=self.__verify, params=post_params, allow_redirects=False, **self._requestArgs)

        if r.status_code == 200 or r.status_code == 302:
            if isJson:
                d = r.json()
                if 'username' in d:
                    del d['username']
                if 'password' in d:
                    del d['password']
                return d
            else:
                return r.text

        # @see self.body for HTTP response body
        raise Exception("Unexpected response from apollo %s: %s" %
                        (r.status_code, r.text))

    def get(self, clientMethod, get_params):
        url = self._wa.apollo_url + self.CLIENT_BASE + clientMethod
        headers = {}

        r = requests.get(url, headers=headers, verify=self.__verify,
                         params=get_params, **self._requestArgs)
        if r.status_code == 200:
            d = r.json()
            if 'username' in d:
                del d['username']
            if 'password' in d:
                del d['password']
            return d
        # @see self.body for HTTP response body
        raise Exception("Unexpected response from apollo %s: %s" %
                        (r.status_code, r.text))


class MetricsClient(Client):
    CLIENT_BASE = '/metrics/'

    def getServerMetrics(self):
        return self.get('metrics', {})


class AnnotationsClient(Client):
    CLIENT_BASE = '/annotationEditor/'

    def _update_data(self, data):
        if not hasattr(self, '_extra_data'):
            raise Exception("Please call setSequence first")

        data.update(self._extra_data)
        return data

    def setSequence(self, sequence, organism):
        self._extra_data = {
            'sequence': sequence,
            'organism': organism,
        }

    def setDescription(self, featureDescriptions):
        data = {
            'features': featureDescriptions,
        }
        data = self._update_data(data)
        return self.request('setDescription', data)

    def setName(self, uniquename, name):
        # TODO
        data = {
            'features': [
                {
                    'uniquename': uniquename,
                    'name': name,
                }
            ],
        }
        data = self._update_data(data)
        return self.request('setName', data)

    def setNames(self, features):
        # TODO
        data = {
            'features': features,
        }
        data = self._update_data(data)
        return self.request('setName', data)

    def setStatus(self, statuses):
        # TODO
        data = {
            'features': statuses,
        }
        data = self._update_data(data)
        return self.request('setStatus', data)

    def setSymbol(self, symbols):
        data = {
            'features': symbols,
        }
        data.update(self._extra_data)
        return self.request('setSymbol', data)

    def getComments(self, features):
        data = {
            'features': features,
        }
        data = self._update_data(data)
        return self.request('getComments', data)

    def addAttribute(self, features):
        data = {
            'features': features,
        }
        data = self._update_data(data)
        return self.request('addAttribute', data)

    def getFeatures(self):
        data = self._update_data({})
        return self.request('getFeatures', data)

    def getSequence(self, uniquename):
        data = {
            'features': [
                {'uniquename': uniquename}
            ]
        }
        data = self._update_data(data)
        return self.request('getSequence', data)

    def addFeature(self, feature, trustme=False):
        if not trustme:
            raise NotImplementedError("Waiting on better docs from project. If you know what you are doing, pass trustme=True to this function.")

        data = {}
        data.update(feature)
        data = self._update_data(data)
        return self.request('addFeature', data)

    def addTranscript(self, transcript, trustme=False):
        if not trustme:
            raise NotImplementedError("Waiting on better docs from project. If you know what you are doing, pass trustme=True to this function.")

        data = {}
        data.update(transcript)
        data = self._update_data(data)
        return self.request('addTranscript', data)

    # addExon, add/delete/updateComments, addTranscript skipped due to docs

    def duplicateTranscript(self, transcriptId):
        data = {
            'features': [{'uniquename': transcriptId}]
        }

        data = self._update_data(data)
        return self.request('duplicateTranscript', data)

    def setTranslationStart(self, uniquename, start):
        data = {
            'features': [{
                'uniquename': uniquename,
                'location': {
                    'fmin': start
                }
            }]
        }
        data = self._update_data(data)
        return self.request('setTranslationStart', data)

    def setTranslationEnd(self, uniquename, end):
        data = {
            'features': [{
                'uniquename': uniquename,
                'location': {
                    'fmax': end
                }
            }]
        }
        data = self._update_data(data)
        return self.request('setTranslationEnd', data)

    def setLongestOrf(self, uniquename):
        data = {
            'features': [{
                'uniquename': uniquename,
            }]
        }
        data = self._update_data(data)
        return self.request('setLongestOrf', data)

    def setBoundaries(self, uniquename, start, end):
        data = {
            'features': [{
                'uniquename': uniquename,
                'location': {
                    'fmin': start,
                    'fmax': end,
                }
            }]
        }
        data = self._update_data(data)
        return self.request('setBoundaries', data)

    def getSequenceAlterations(self):
        data = {
        }
        data = self._update_data(data)
        return self.request('getSequenceAlterations', data)

    def setReadthroughStopCodon(self, uniquename):
        data = {
            'features': [{
                'uniquename': uniquename,
            }]
        }
        data = self._update_data(data)
        return self.request('setReadthroughStopCodon', data)

    def deleteSequenceAlteration(self, uniquename):
        data = {
            'features': [{
                'uniquename': uniquename,
            }]
        }
        data = self._update_data(data)
        return self.request('deleteSequenceAlteration', data)

    def flipStrand(self, uniquenames):
        data = {
            'features': [
                {'uniquename': x} for x in uniquenames
            ]
        }
        data = self._update_data(data)
        return self.request('flipStrand', data)

    def mergeExons(self, exonA, exonB):
        data = {
            'features': [
                {'uniquename': exonA},
                {'uniquename': exonB},
            ]
        }
        data = self._update_data(data)
        return self.request('mergeExons', data)

    # def splitExon(): pass

    def deleteFeatures(self, uniquenames):
        assert isinstance(uniquenames, collections.Iterable)
        data = {
            'features': [
                {'uniquename': x} for x in uniquenames
            ]
        }
        data = self._update_data(data)
        return self.request('deleteFeature', data)

    # def deleteExon(): pass

    # def makeIntron(self, uniquename, ): pass

    def getSequenceSearchTools(self):
        return self.get('getSequenceSearchTools', {})

    def getCannedComments(self):
        return self.get('getCannedComments', {})

    def searchSequence(self, searchTool, sequence, database):
        data = {
            'key': searchTool,
            'residues': sequence,
            'database_id': database,
        }
        return self.request('searchSequences', data)

    def getGff3(self, uniquenames):
        assert isinstance(uniquenames, collections.Iterable)
        data = {
            'features': [
                {'uniquename': x} for x in uniquenames
            ]
        }
        data = self._update_data(data)
        return self.request('getGff3', data, isJson=False)


class GroupsClient(Client):
    CLIENT_BASE = '/group/'

    def createGroup(self, name):
        data = {'name': name}
        return self.request('createGroup', data)

    def getOrganismPermissionsForGroup(self, group):
        data = {
            'id': group.groupId,
            'name': group.name,
        }
        return self.request('getOrganismPermissionsForGroup', data)

    def loadGroups(self, group=None):
        data = {}
        if group is not None:
            data['groupId'] = group.groupId

        return self.request('loadGroups', data)

    def deleteGroup(self, group):
        data = {
            'id': group.groupId,
            'name': group.name,
        }
        return self.request('deleteGroup', data)

    def updateGroup(self, group, newName):
        # TODO: Sure would be nice if modifying ``group.name`` would invoke
        # this?
        data = {
            'id': group.groupId,
            'name': newName,
        }
        return self.request('updateGroup', data)

    def updateOrganismPermission(self, group, organismName,
                                 administrate=False, write=False, read=False,
                                 export=False):
        data = {
            'groupId': group.groupId,
            'name': organismName,
            'administrate': administrate,
            'write': write,
            'export': export,
            'read': read,
        }
        return self.request('updateOrganismPermission', data)

    def updateMembership(self, group, users):
        data = {
            'groupId': group.groupId,
            'user': [user.email for user in users]
        }
        return self.request('updateMembership', data)


class IOClient(Client):
    CLIENT_BASE = '/IOService/'

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

        return self.request('write', data, isJson=output == 'file')

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

    def addOrganism(self, commonName, directory, blatdb=None, species=None,
                    genus=None, public=False):
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

    def findOrganismByCn(self, cn):
        orgs = self.findAllOrganisms()
        orgs = [x for x in orgs if x['commonName'] == cn]
        if len(orgs) == 0:
            raise Exception("Unknown common name")
        else:
            return orgs[0]

    def findOrganismById(self, id_number):
        orgs = self.findAllOrganisms()
        orgs = [x for x in orgs if str(x['id']) == str(id_number)]
        if len(orgs) == 0:
            raise Exception("Unknown ID")
        else:
            return orgs[0]

    def deleteOrganism(self, organismId):
        return self.request('deleteOrganism', {'id': organismId})

    def deleteOrganismFeatures(self, organismId):
        return self.request('deleteOrganismFeatures', {'id': organismId})

    def getSequencesForOrganism(self, commonName):
        return self.request('getSequencesForOrganism', {'organism': commonName})

    def updateOrganismInfo(self, organismId, commonName, directory, blatdb=None, species=None, genus=None, public=False):
        data = {
            'id': organismId,
            'name': commonName,
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

    # Real one
    # def getOrganismPermissionsForUser(self, user):
    # data = {
    # 'userId': user.userId,
    # }
    # return self.request('getOrganismPermissionsForUser', data)

    # Utter frigging hack
    def getOrganismPermissionsForUser(self, user):
        return self.loadUser(user).organismPermissions

    def updateOrganismPermission(self, user, organism, administrate=False,
                                 write=False, export=False, read=False):
        data = {
            'userId': user.userId,
            'organism': organism,
            'ADMINISTRATE': administrate,
            'WRITE': write,
            'EXPORT': export,
            'READ': read,
        }
        return self.request('updateOrganismPermission', data)

    def loadUser(self, user):
        return self.loadUserById(user.userId)

    def loadUserById(self, userId):
        res = self.request('loadUsers', {'userId': userId})
        if isinstance(res, list):
            # We can only match one, right?
            return UserObj(**res[0])
        else:
            return res

    def loadUsers(self, email=None):
        res = self.request('loadUsers', {})
        data = [UserObj(**x) for x in res]
        if email is not None:
            data = [x for x in data if x.username == email]

        return data

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


class RemoteRecord(Client):
    CLIENT_BASE = None

    def ParseRecord(self, cn):
        org = self._wa.organisms.findOrganismByCn(cn)
        self._wa.annotations.setSequence(org['commonName'], org['id'])

        data = StringIO.StringIO(self._wa.io.write(
            exportType='GFF3',
            seqType='genomic',
            exportAllSequences=False,
            exportGff3Fasta=True,
            output="text",
            exportFormat="text",
            sequences=cn,
        ))
        data.seek(0)

        for record in GFF.parse(data):
            yield WebApolloSeqRecord(record, self._wa)


class WebApolloSeqRecord(object):

    def __init__(self, sr, wa):
        self._sr = sr
        self._wa = wa

    def __dir__(self):
        return dir(self._sr)

    def __getattr__(self, key):
        if key in ('_sr', '_wa'):
            return self.__dict__[key]
        else:
            if key == 'features':
                return (WebApolloSeqFeature(x, self._wa)
                        for x in self._sr.__dict__[key])
            else:
                return self._sr.__dict__[key]

    def __setattr__(self, key, value):
        if key in ('_sd', '_wa'):
            self.__dict__[key] = value
        else:
            self._sr.__dict__[key] = value
            # Methods acting on the SeqRecord object


class WebApolloSeqFeature(object):

    def __init__(self, sf, wa):
        self._sf = sf
        self._wa = wa

    def __dir__(self):
        return dir(self._sf)

    def __getattr__(self, key):
        if key in ('_sf', '_wa'):
            return self.__dict__[key]
        else:
            return self._sf.__dict__[key]

    def __setattr__(self, key, value):
        if key in ('_sf', '_wa'):
            self.__dict__[key] = value
        else:
            # Methods acting on the SeqFeature object
            if key == 'location':
                if value.strand != self._sf.location.strand:
                    self.wa.annotations.flipStrand(
                        self._sf.qualifiers['ID'][0]
                    )

                self.wa.annotations.setBoundaries(
                    self._sf.qualifiers['ID'][0],
                    value.start,
                    value.end,
                )

                self._sf.__dict__[key] = value
            else:
                self._sf.__dict__[key] = value


def _tnType(feature):
    if feature.type in ('gene', 'mRNA', 'exon', 'CDS'):
        return feature.type
    else:
        return 'exon'


def _yieldFeatData(features):
    for f in features:
        current = {
            'location': {
                'strand': f.strand,
                'fmin': int(f.location.start),
                'fmax': int(f.location.end),
            },
            'type': {
                'name': _tnType(f),
                'cv': {
                    'name': 'sequence',
                }
            },
        }
        if f.type in ('gene', 'mRNA'):
            current['name'] = f.qualifiers.get('Name', [f.id])[0]
        if hasattr(f, 'sub_features') and len(f.sub_features) > 0:
            current['children'] = [x for x in _yieldFeatData(f.sub_features)]

        yield current


def featuresToFeatureSchema(features):
    compiled = []
    for feature in features:
        if feature.type != 'gene':
            log.warn("Not able to handle %s features just yet...", feature.type)
            continue

        for x in _yieldFeatData([feature]):
            compiled.append(x)
    return compiled


def accessible_organisms(user, orgs):
    permissionMap = {
        x['organism']: x['permissions']
        for x in user.organismPermissions
        if user.isAdmin() or any(PERM in x['permissions'] for PERM in ('WRITE', 'READ', 'ADMINISTRATE'))
    }

    for org in sorted(orgs, key=lambda x: x.get('commonName', '')):
        if org.get('commonName', None) in permissionMap:
            yield org


def galaxy_list_orgs(trans, *args, **kwargs):
    email = trans.get_user().email

    wa = WebApolloInstance(
        os.environ.get('GALAXY_WEBAPOLLO_URL', 'https://example.com'),
        os.environ.get('GALAXY_WEBAPOLLO_USER', 'admin'),
        os.environ.get('GALAXY_WEBAPOLLO_PASSWORD', 'admin')
    )

    # Get the user
    gx_user = AssertUser(wa.users.loadUsers(email=email))
    # Find ALL organisms
    all_orgs = wa.organisms.findAllOrganisms()
    # And then filter by those which the user has permissions on.
    orgs = accessible_organisms(gx_user, all_orgs)

    response = []
    for org in orgs:
        if 'genus' in org and org['genus']:
            response.append(("{commonName} ({genus} {species})".format(**org), org['id'], False))
        else:
            response.append(("{commonName}".format(**org), org['id'], False))
    print(response)
    return response
