#!/usr/bin/env python
from galaxygetopt.ggo import GalaxyGetOpt as GGO
import sys
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()
from Bio import Entrez
from Bio import SeqIO
import tempfile
import os


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
            [
                'data',
                'Exported data',
                {
                    'validate': 'File/Output',
                    'required': True,
                    'default': 'export',
                    'data_format': 'genomic/annotated',
                    'default_format': 'Genbank',
                }
            ]
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
    import re
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
                id = record.id
                out_name = os.path.join("gbk_out", '%s.gbk' % id)

                # Collision free names
                if os.path.exists(out_name):
                    i = 0
                    while os.path.exists(os.path.join("gbk_out", '%s.x%s.gbk' % (id, i))):
                        i += 1
                    out_name = os.path.join("gbk_out", '%s.x%s.gbk' % (id, i))

                with open(out_name, 'w') as handle:
                    SeqIO.write(record, handle, 'genbank')
    else:
        for record_set in results:
            out_name = os.path.join("gbk_out", 'merged.gbk')
            with open(out_name, 'a') as handle:
                for record in record_set:
                    SeqIO.write(record, handle, 'genbank')
