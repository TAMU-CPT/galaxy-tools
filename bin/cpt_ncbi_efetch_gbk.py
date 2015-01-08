#!/usr/bin/env python
from galaxygetopt.ggo import GalaxyGetOpt as GGO
import sys
import logging
import string
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()
from Bio import Entrez
from Bio import SeqIO
import tempfile
import os
import re
phage_in_middle = re.compile('^(?P<host>.*)\s*phage (?P<phage>.*)$')
bacteriophage_in_middle = re.compile('^(?P<host>.*)\s*bacteriophage (?P<phage>.*)$')
starts_with_phage = re.compile('^(bacterio|vibrio|Bacterio|Vibrio|)?[Pp]hage (?P<phage>.*)$')
new_style_names = re.compile('(?P<phage>v[A-Z]_[A-Z][a-z]{2}_.*)')

valid_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)


__doc__ = """
NCBI Entrez Fetcher
=======================

fetch records from NCBI via Entrez EFetch
"""


def efetch(email=None, accessions=None, db="nucleotide", query_pagination_size=500):
    if accessions is None:
        log.error("No accessions provided")
        sys.exit(1)
    if not isinstance(accessions, list):
        log.error("Accessions must be a list")
        sys.exit(1)
    if len(accessions) == 0:
        log.error("No accessions provided")
        sys.exit(1)
    if accessions is None:
        log.error("No email provided")
        sys.exit(1)

    Entrez.email = email
    Entrez.tool = "Galaxy"
    cleanup = []

    records = []

    for sublist in [accessions[i:i+query_pagination_size] for i in range(0, len(accessions), query_pagination_size)]:
        handle = Entrez.efetch(db=db, id=sublist, rettype="gb",
                               retmode="text")
        # Serialize
        tmp = tempfile.NamedTemporaryFile(delete=False)
        tmp.write(handle.read())
        tmp.close()
        # SeqIO.read works with the handle, however it doesn't accept multiple
        # records, so we have to use SeqIO.parse which seems to want a file
        # handle
        record = SeqIO.parse(tmp.name, "genbank")
        handle.close()
        cleanup.append(tmp.name)
        records.append(record)

    return records


def db_finder():
    Entrez.email = "esr@tamu.edu"
    Entrez.tool = "Galaxy"
    data = {}

    handle = Entrez.einfo()
    result = Entrez.read(handle)
    handle.close()

    for db in result['DbList']:
    #for db in ['unigene']:
        data[db] = {}
        try:
            h = Entrez.einfo(db=db)
            r = Entrez.read(h)
            h.close()
            data[db] = r
        except:
            pass

    import yaml
    with open('out.yml', 'w') as handle:
        handle.write(yaml.dump(data))


def name_parser(name):
    host = None
    phage = None
    name = name.replace(', complete genome', '')

    m = bacteriophage_in_middle.match(name)
    if m:
        host = m.group('host')
        phage = m.group('phage')
        return (host, phage)

    m = phage_in_middle.match(name)
    if m:
        host = m.group('host')
        phage = m.group('phage')
        return (host, phage)

    m = starts_with_phage.match(name)
    if m:
        phage = m.group('phage')
        return (host, phage)

    m = new_style_names.match(name)
    if m:
        phage = m.group('phage')
        return (host, phage)

    return (host, phage)

if __name__ == '__main__':
    opts = GGO(
        options=[
            ['accessions_list', 'Comma/space/tab separated accessions', {'validate': 'String'}],
            ['accessions_file', 'Comma/space/tab separated accessions (from File)', {'validate': 'File/Input'}],
            ['email', 'User email address', {'validate': 'String'}],
            ['split', 'Split genbank record (for multiple records)', {'validate': 'Flag'}],
            ['db', 'Database to fetch from', {'validate': 'String', 'default': 'nucleotide'}],
        ],
        outputs=[
        ],
        defaults={
            'appid': 'edu.tamu.cpt.ncbi.entrez.efetch',
            'appname': 'NCBI Entrez EFetch',
            'appvers': '1.0.0',
            'appdesc': 'downloads data from NCBI via Entrez EFetch',
        },
        tests=[],
        doc=__doc__
    )
    options = opts.params()

    al = options['accessions_list']
    af = options['accessions_file']

    accessions = []
    if al is not None and al is not "None":
        for value in re.findall(r'[A-Za-z0-9_-]+', al):
            accessions.append(value)

    if af is not None and af is not "None":
        data = af.read()
        for value in re.findall(r'[A-Za-z0-9_-]+', data):
            accessions.append(value)

    results = efetch(email=options['email'], accessions=accessions, db=options['db'])

    try:
        os.mkdir("gbk_out")
    except:
        pass

    if options['split']:
        for record_set in results:
            for record in record_set:
                (host, phage) = name_parser(record.description)
                if host is None:
                    host = ""
                if phage is None:
                    phage = ""
                # Move to regex
                host = host.strip()
                phage = phage.strip()

                if host and phage:
                    id = "%s [%s]" % (phage, host)
                elif not host and phage:
                    id = "%s" % (phage)
                else:
                    id = record.id

                id_fixed = ''.join(c for c in id if c in valid_chars)

                out_name = os.path.join("gbk_out", '%s.gbk' % id_fixed)

                # Collision free names
                if os.path.exists(out_name):
                    i = 0
                    while os.path.exists(os.path.join("gbk_out", '%s.x%s.gbk' % (id_fixed, i))):
                        i += 1
                    out_name = os.path.join("gbk_out", '%s.x%s.gbk' % (id_fixed, i))

                with open(out_name, 'w') as handle:
                    SeqIO.write(record, handle, 'genbank')


    else:
        for record_set in results:
            out_name = os.path.join("gbk_out", 'merged.gbk')
            with open(out_name, 'a') as handle:
                for record in record_set:
                    SeqIO.write(record, handle, 'genbank')
