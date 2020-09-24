from __future__ import print_function

import argparse
import collections
import json
import logging
import os
import time
from abc import abstractmethod

import requests

from six.moves.builtins import next
from six.moves.builtins import object

import yaml

logging.getLogger("requests").setLevel(logging.CRITICAL)
log = logging.getLogger()


#############################################
#      BEGIN IMPORT OF CACHING LIBRARY      #
#############################################
# This code is licensed under the MIT       #
# License and is a copy of code publicly    #
# available in rev.                         #
# e27332bc82f4e327aedaec17c9b656ae719322ed  #
# of https://github.com/tkem/cachetools/    #
#############################################

class DefaultMapping(collections.MutableMapping):

    __slots__ = ()

    @abstractmethod
    def __contains__(self, key):  # pragma: nocover
        return False

    @abstractmethod
    def __getitem__(self, key):  # pragma: nocover
        if hasattr(self.__class__, '__missing__'):
            return self.__class__.__missing__(self, key)
        else:
            raise KeyError(key)

    def get(self, key, default=None):
        if key in self:
            return self[key]
        else:
            return default

    __marker = object()

    def pop(self, key, default=__marker):
        if key in self:
            value = self[key]
            del self[key]
        elif default is self.__marker:
            raise KeyError(key)
        else:
            value = default
        return value

    def setdefault(self, key, default=None):
        if key in self:
            value = self[key]
        else:
            self[key] = value = default
        return value


DefaultMapping.register(dict)


class _DefaultSize(object):
    def __getitem__(self, _):
        return 1

    def __setitem__(self, _, value):
        assert value == 1

    def pop(self, _):
        return 1


class Cache(DefaultMapping):
    """Mutable mapping to serve as a simple cache or cache base class."""

    __size = _DefaultSize()

    def __init__(self, maxsize, missing=None, getsizeof=None):
        if missing:
            self.__missing = missing
        if getsizeof:
            self.__getsizeof = getsizeof
            self.__size = dict()
        self.__data = dict()
        self.__currsize = 0
        self.__maxsize = maxsize

    def __repr__(self):
        return '%s(%r, maxsize=%r, currsize=%r)' % (
            self.__class__.__name__,
            list(self.__data.items()),
            self.__maxsize,
            self.__currsize,
        )

    def __getitem__(self, key):
        try:
            return self.__data[key]
        except KeyError:
            return self.__missing__(key)

    def __setitem__(self, key, value):
        maxsize = self.__maxsize
        size = self.getsizeof(value)
        if size > maxsize:
            raise ValueError('value too large')
        if key not in self.__data or self.__size[key] < size:
            while self.__currsize + size > maxsize:
                self.popitem()
        if key in self.__data:
            diffsize = size - self.__size[key]
        else:
            diffsize = size
        self.__data[key] = value
        self.__size[key] = size
        self.__currsize += diffsize

    def __delitem__(self, key):
        size = self.__size.pop(key)
        del self.__data[key]
        self.__currsize -= size

    def __contains__(self, key):
        return key in self.__data

    def __missing__(self, key):
        value = self.__missing(key)
        try:
            self.__setitem__(key, value)
        except ValueError:
            pass  # value too large
        return value

    def __iter__(self):
        return iter(self.__data)

    def __len__(self):
        return len(self.__data)

    @staticmethod
    def __getsizeof(value):
        return 1

    @staticmethod
    def __missing(key):
        raise KeyError(key)

    @property
    def maxsize(self):
        """The maximum size of the cache."""
        return self.__maxsize

    @property
    def currsize(self):
        """The current size of the cache."""
        return self.__currsize

    def getsizeof(self, value):
        """Return the size of a cache element's value."""
        return self.__getsizeof(value)


class _Link(object):

    __slots__ = ('key', 'expire', 'next', 'prev')

    def __init__(self, key=None, expire=None):
        self.key = key
        self.expire = expire

    def __reduce__(self):
        return _Link, (self.key, self.expire)

    def unlink(self):
        next = self.next
        prev = self.prev
        prev.next = next
        next.prev = prev


class _Timer(object):

    def __init__(self, timer):
        self.__timer = timer
        self.__nesting = 0

    def __call__(self):
        if self.__nesting == 0:
            return self.__timer()
        else:
            return self.__time

    def __enter__(self):
        if self.__nesting == 0:
            self.__time = time = self.__timer()
        else:
            time = self.__time
        self.__nesting += 1
        return time

    def __exit__(self, *exc):
        self.__nesting -= 1

    def __reduce__(self):
        return _Timer, (self.__timer,)

    def __getattr__(self, name):
        return getattr(self.__timer, name)


class TTLCache(Cache):
    """LRU Cache implementation with per-item time-to-live (TTL) value."""

    def __init__(self, maxsize, ttl, timer=time.time, missing=None,
                 getsizeof=None):
        Cache.__init__(self, maxsize, missing, getsizeof)
        self.__root = root = _Link()
        root.prev = root.next = root
        self.__links = collections.OrderedDict()
        self.__timer = _Timer(timer)
        self.__ttl = ttl

    def __contains__(self, key):
        try:
            link = self.__links[key]  # no reordering
        except KeyError:
            return False
        else:
            return not (link.expire < self.__timer())

    def __getitem__(self, key, cache_getitem=Cache.__getitem__):
        try:
            link = self.__getlink(key)
        except KeyError:
            expired = False
        else:
            expired = link.expire < self.__timer()
        if expired:
            return self.__missing__(key)
        else:
            return cache_getitem(self, key)

    def __setitem__(self, key, value, cache_setitem=Cache.__setitem__):
        with self.__timer as time:
            self.expire(time)
            cache_setitem(self, key, value)
        try:
            link = self.__getlink(key)
        except KeyError:
            self.__links[key] = link = _Link(key)
        else:
            link.unlink()
        link.expire = time + self.__ttl
        link.next = root = self.__root
        link.prev = prev = root.prev
        prev.next = root.prev = link

    def __delitem__(self, key, cache_delitem=Cache.__delitem__):
        cache_delitem(self, key)
        link = self.__links.pop(key)
        link.unlink()
        if link.expire < self.__timer():
            raise KeyError(key)

    def __iter__(self):
        root = self.__root
        curr = root.next
        while curr is not root:
            # "freeze" time for iterator access
            with self.__timer as time:
                if not (curr.expire < time):
                    yield curr.key
            curr = curr.next

    def __len__(self):
        root = self.__root
        curr = root.next
        time = self.__timer()
        count = len(self.__links)
        while curr is not root and curr.expire < time:
            count -= 1
            curr = curr.next
        return count

    def __setstate__(self, state):
        self.__dict__.update(state)
        root = self.__root
        root.prev = root.next = root
        for link in sorted(self.__links.values(), key=lambda obj: obj.expire):
            link.next = root
            link.prev = prev = root.prev
            prev.next = root.prev = link
        self.expire(self.__timer())

    def __repr__(self, cache_repr=Cache.__repr__):
        with self.__timer as time:
            self.expire(time)
            return cache_repr(self)

    @property
    def currsize(self):
        with self.__timer as time:
            self.expire(time)
            return super(TTLCache, self).currsize

    @property
    def timer(self):
        """The timer function used by the cache."""
        return self.__timer

    @property
    def ttl(self):
        """The time-to-live value of the cache's items."""
        return self.__ttl

    def expire(self, time=None):
        """Remove expired items from the cache."""
        if time is None:
            time = self.__timer()
        root = self.__root
        curr = root.next
        links = self.__links
        cache_delitem = Cache.__delitem__
        while curr is not root and curr.expire < time:
            cache_delitem(self, curr.key)
            del links[curr.key]
            next = curr.next
            curr.unlink()
            curr = next

    def clear(self):
        with self.__timer as time:
            self.expire(time)
            Cache.clear(self)

    def get(self, *args, **kwargs):
        with self.__timer:
            return Cache.get(self, *args, **kwargs)

    def pop(self, *args, **kwargs):
        with self.__timer:
            return Cache.pop(self, *args, **kwargs)

    def setdefault(self, *args, **kwargs):
        with self.__timer:
            return Cache.setdefault(self, *args, **kwargs)

    def popitem(self):
        """Remove and return the `(key, value)` pair least recently used that
        has not already expired.

        """
        with self.__timer as time:
            self.expire(time)
            try:
                key = next(iter(self.__links))
            except StopIteration:
                raise KeyError('%s is empty' % self.__class__.__name__)
            else:
                return (key, self.pop(key))

    if hasattr(collections.OrderedDict, 'move_to_end'):
        def __getlink(self, key):
            value = self.__links[key]
            self.__links.move_to_end(key)
            return value
    else:
        def __getlink(self, key):
            value = self.__links.pop(key)
            self.__links[key] = value
            return value


#############################################
#       END IMPORT OF CACHING LIBRARY       #
#############################################


cache = TTLCache(
    100,  # Up to 100 items
    5 * 60  # 5 minute cache life
)
userCache = TTLCache(
    2,  # Up to 2 items
    60  # 1 minute cache life
)


class UnknownUserException(Exception):
    pass


def WAAuth(parser):
    parser.add_argument('apollo', help='Complete Apollo URL')
    parser.add_argument('username', help='WA Username')
    parser.add_argument('password', help='WA Password')


def OrgOrGuess(parser):
    parser.add_argument('--org_json', type=argparse.FileType("r"), help='Apollo JSON output, source for common name')
    parser.add_argument('--org_raw', help='Common Name')
    parser.add_argument('--org_id', help='Organism ID')


def CnOrGuess(parser):
    OrgOrGuess(parser)
    parser.add_argument('--seq_fasta', type=argparse.FileType("r"), help='Fasta file, IDs used as sequence sources')
    parser.add_argument('--seq_raw', nargs='*', help='Sequence Names')


def AssertUser(user_list):
    if len(user_list) == 0:
        raise UnknownUserException()
    elif len(user_list) == 1:
        return user_list[0]
    else:
        raise Exception("Too many users!")


class WebApolloInstance(object):

    def __init__(self):

        if 'ARROW_GLOBAL_CONFIG_PATH' in os.environ:

            with open(os.environ['ARROW_GLOBAL_CONFIG_PATH'], 'r') as config:
                conf = yaml.safe_load(config)
                try:
                    instance_name = conf['__default']
                except KeyError:
                    raise Exception("Unknown Apollo instance and no __default provided")
                self.apollo_url = conf[instance_name]['url']
                self.username = conf[instance_name]['username']
                self.password = conf[instance_name]['password']
        else:
            self.apollo_url = os.environ['GALAXY_WEBAPOLLO_URL']
            self.username = os.environ['GALAXY_WEBAPOLLO_USER']
            self.password = os.environ['GALAXY_WEBAPOLLO_PASSWORD']

        self.groups = GroupsClient(self)
        self.organisms = OrganismsClient(self)
        self.users = UsersClient(self)

    def __str__(self):
        return '<WebApolloInstance at %s>' % self.apollo_url

    def requireUser(self, email):
        cacheKey = 'user-list'
        try:
            # Get the cached value
            data = userCache[cacheKey]
        except KeyError:
            # If we hit a key error above, indicating that
            # we couldn't find the key, we'll simply re-request
            # the data
            data = self.users.loadUsers()
            userCache[cacheKey] = data

        return AssertUser([x for x in data if x.username == email])


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

        r = requests.post(url, data=json.dumps(data), headers=headers,
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


class GroupsClient(Client):
    CLIENT_BASE = '/group/'

    def loadGroups(self, group=None):
        res = self.request('loadGroups', {})
        data = [GroupObj(**x) for x in res]
        if group is not None:
            data = [x for x in data if x.name == group]

        return data


class OrganismsClient(Client):
    CLIENT_BASE = '/organism/'

    def findAllOrganisms(self):
        orgs = self.request('findAllOrganisms', {})
        if not isinstance(orgs, (list,)):
            orgs = []
        return orgs


class UsersClient(Client):
    CLIENT_BASE = '/user/'

    def loadUsers(self):
        res = self.request('loadUsers', {})

        data = [UserObj(**x) for x in res]

        return data


def handle_credentials(user):
    if hasattr(user, 'new_password'):
        f = open("Apollo_credentials.txt", "w")
        f.write('Username:\t%s\nPassword:\t%s' % (user.username, user.new_password))


def accessible_organisms(user, orgs):
    permissionMap = {
        x['organism']: x['permissions']
        for x in user.organismPermissions
        if 'WRITE' in x['permissions'] or 'READ' in x['permissions'] or 'ADMINISTRATE' in x['permissions'] or user.role == 'ADMIN'
    }

    if 'error' in orgs:
        raise Exception("Error received from Apollo server: \"%s\"" % orgs['error'])

    return [
        (org['commonName'], org['id'], False)
        for org in orgs#sorted(orgs, key=lambda x: x['commonName'])
        if org['commonName'] in permissionMap
    ]


def galaxy_list_groups(trans, *args, **kwargs):
    email = trans.get_user().email
    wa = WebApolloInstance()

    # Key for cached data
    cacheKey = 'groups-' + email
    # We don't want to trust "if key in cache" because between asking and fetch
    # it might through key error.
    if cacheKey not in cache:
        # However if it ISN'T there, we know we're safe to fetch + put in
        # there.
        data = _galaxy_list_groups(wa, *args, **kwargs)
        cache[cacheKey] = data
        return data
    try:
        # The cache key may or may not be in the cache at this point, it
        # /likely/ is. However we take no chances that it wasn't evicted between
        # when we checked above and now, so we reference the object from the
        # cache in preparation to return.
        data = cache[cacheKey]
        return data
    except KeyError:
        # If access fails due to eviction, we will fail over and can ensure that
        # data is inserted.
        data = _galaxy_list_groups(wa, *args, **kwargs)
        cache[cacheKey] = data
        return data


def _galaxy_list_groups(wa, *args, **kwargs):
    # Fetch the groups.
    group_data = []
    for group in wa.groups.loadGroups():
        # Reformat
        group_data.append((group.name, group.name, False))
    return group_data


def galaxy_list_users(trans, *args, **kwargs):
   # email = trans.get_user().email
    wa = WebApolloInstance()
    userDicts = wa.users.loadUsers()
    outList = []
    #return "String [0]: " + str(userDicts[0]) + "\n------\nDir: " + str(dir(userDicts[0]))
    for x in userDicts:
        outList.append((x.username, x.username, False))
    return outList

def galaxy_list_orgs(trans, *args, **kwargs):
    email = trans.get_user().email
    wa = WebApolloInstance()
    try:
        gx_user = wa.requireUser(email)
    except UnknownUserException:
        return []

    # Key for cached data
    cacheKey = 'orgs-' + email
    if cacheKey not in cache:
        data = _galaxy_list_orgs(wa, gx_user, *args, **kwargs)
        cache[cacheKey] = data
        return data
    try:
        data = cache[cacheKey]
        return data
    except KeyError:
        data = _galaxy_list_orgs(wa, gx_user, *args, **kwargs)
        cache[cacheKey] = data
        return data


def _galaxy_list_orgs(wa, gx_user, *args, **kwargs):
    # Fetch all organisms
    all_orgs = wa.organisms.findAllOrganisms()
    # Figure out which are accessible to the user
    orgs = accessible_organisms(gx_user, all_orgs)
    # Return org list
    return orgs


# This is all for implementing the command line interface for testing.
class obj(object):
    pass


class fakeTrans(object):

    def __init__(self, username):
        self.un = username

    def get_user(self):
        o = obj()
        o.email = self.un
        return o


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Test access to apollo server')
    parser.add_argument('email', help='Email of user to test')
    parser.add_argument('--action', choices=['org', 'group'], default='org', help='Data set to test, fetch a list of groups or orgs known to the requesting user.')
    args = parser.parse_args()

    trans = fakeTrans(args.email)
    if args.action == 'org':
        print(galaxy_list_orgs(trans))
    elif args.action == 'group':
        print(galaxy_list_groups(trans))
