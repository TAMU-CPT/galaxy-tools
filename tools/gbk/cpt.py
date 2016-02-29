import re

PHAGE_IN_MIDDLE = re.compile('^(?P<host>.*)\s*phage (?P<phage>.*)$')
BACTERIOPHAGE_IN_MIDDLE = re.compile('^(?P<host>.*)\s*bacteriophage (?P<phage>.*)$')
STARTS_WITH_PHAGE = re.compile('^(bacterio|vibrio|Bacterio|Vibrio|)?[Pp]hage (?P<phage>.*)$')
NEW_STYLE_NAMES = re.compile('(?P<phage>v[A-Z]_[A-Z][a-z]{2}_.*)')

def phage_name_parser(name):
    host = None
    phage = None
    name = name.replace(', complete genome.', '')
    name = name.replace(', complete genome', '')

    m = BACTERIOPHAGE_IN_MIDDLE.match(name)
    if m:
        host = m.group('host')
        phage = m.group('phage')
        return (host, phage)

    m = PHAGE_IN_MIDDLE.match(name)
    if m:
        host = m.group('host')
        phage = m.group('phage')
        return (host, phage)

    m = STARTS_WITH_PHAGE.match(name)
    if m:
        phage = m.group('phage')
        return (host, phage)

    m = NEW_STYLE_NAMES.match(name)
    if m:
        phage = m.group('phage')
        return (host, phage)

    return (host, phage)


