#!/usr/bin/env python
import subprocess
import argparse
import os
import yaml
import re
from Bio import SeqIO


SBT_TPL = '''Seq-submit ::= {
  contact {
    contact {
      name name {
%s
      },
      affil std {
%s
      }
    }
  },
  cit {
    authors {
      names std {
%s
      },
      affil std {
%s
      }
    }
  },
  subtype new
}'''

SBT_TPL2 = '''
Seqdesc ::= pub {
  pub {
    gen {
      cit "%s",
      authors {
        names std {
          {
%s
          }
        },
        affil std {
%s
        }
      },
      title "%s"
    }
  }
}
'''


def create_contact_string(offset, author_dn_list, people):
    output = []
    for dn in author_dn_list:
        person = people['people'][dn]
        letters = [person['givenName'][0]]
        if 'initials' in person:
            unparsed_suffix = re.sub('[^A-Za-z]', '', person['initials'])
            letters += list(unparsed_suffix)

        data = {
            'last': person['sn'],
            'first': person['givenName'],
            'initials': '.'.join([x.upper() for x in letters]),
            'suffix': person.get('generationQualifier', ''),
        }
        data_str = [
            '  ' * offset + '{',
            '  ' * (offset + 1) + 'name name {',
            create_string(offset + 2, data, ['last', 'first', 'initials', 'suffix']),
            '  ' * (offset + 1) + '}',
            '  ' * offset + '}'
        ]
        output.append('\n'.join(data_str))
    return ',\n'.join(output)


def create_affil_string(offset, dn, people):
    org = people['orgs'][dn]
    affil_info = {
        'affil': org['ou'],
        'div': '',
        'city': org['l'],
        'sub': org['st'],
        'country': 'USA',
        'street': org['street'],
        'email': 'cpt@tamu.edu',
        'fax': '',
        'phone': org['telephoneNumber'],
        'zip': org['postalCode']
    }
    afk = ['affil', 'div', 'city', 'sub', 'country', 'street', 'email', 'fax', 'phone', 'zip']
    data = [
        '  ' * offset + '{',
        '  ' * (offset + 1) + 'name name {',
        create_string(offset + 2, affil_info, afk),
        '  ' * (offset + 1) + '}',
        '  ' * offset + '}',

    ]
    return '\n'.join(data)


def create_string(length, data, keys):
    data = [
        ('  ' * length) + '%s "%s"' % (key, data[key])
        for key in keys
    ]
    return ",\n".join(data)
    # return join(",\n",map { '  ' x $length . $_ . ' "' . ($hash{$_} ) . '"' } @arr);


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--genbank_submission_title')
    parser.add_argument('--genbank_record_author', nargs='+')
    parser.add_argument('--genbank_record_contact', nargs='+')
    parser.add_argument('--genbank_record_affiliation')
    parser.add_argument('--paper_submission_title')
    parser.add_argument('--paper_record_author', nargs='+')
    parser.add_argument('--paper_record_contact', nargs='+')
    parser.add_argument('--paper_record_affiliation')
    parser.add_argument('--genbank_file', type=argparse.FileType("r"))
    parser.add_argument('--people', type=argparse.FileType("r"))
    parser.add_argument('--tmpdir')

    args = parser.parse_args()

    os.makedirs(args.tmpdir)
    tmpdir = args.tmpdir

    directory = yaml.load(args.people)
    record = SeqIO.read(args.genbank_file, 'genbank')
    basename = os.path.join(tmpdir, args.genbank_submission_title)

    # Fasta file
    with open(basename + '.fsa', 'w') as handle:
        handle.write(">{0} [organism=Phage {0}]\n".format(args.genbank_submission_title))
        seq = str(record.seq)
        for i in range(0, len(seq), 80):
            handle.write(seq[i:i + 80] + '\n')

    # Feature table
    with open(basename + '.tbl', 'w') as handle:
        handle.write('>Feature %s\n' % args.genbank_submission_title)
        for feature in record.features:
            if feature.type == 'source':
                continue

            if feature.strand == -1:
                handle.write("{0}\t{1}\t{2}\n".format(
                    feature.location.end,
                    feature.location.start + 1,
                    feature.type
                ))
            else:
                handle.write("{0}\t{1}\t{2}\n".format(
                    feature.location.start + 1,
                    feature.location.end,
                    feature.type
                ))

            for qualifier in feature.qualifiers:
                for value in feature.qualifiers[qualifier]:
                    if qualifier == 'source':
                        continue
                    handle.write('\t\t\t%s\t%s\n' % (qualifier, value))

    # SBT
    with open(basename + '.sbt', 'w') as handle:
        contact_string = create_contact_string(4, args.genbank_record_contact, directory)
        affil_string = create_affil_string(4, args.genbank_record_affiliation, directory)
        affil_string2 = create_affil_string(5, args.paper_record_affiliation, directory)
        record_author_str = create_contact_string(4, args.genbank_record_author, directory)
        paper_author_str = create_contact_string(6, args.paper_record_author, directory)

        publication_status = 'unpublished'
        publication_title = args.paper_submission_title

        handle.write(SBT_TPL % (contact_string, affil_string, record_author_str, affil_string))
        handle.write(SBT_TPL2 % (publication_status, paper_author_str, affil_string2, publication_title))

    # tbl2asn
    subprocess.check_call([
        'tbl2asn', '-p', tmpdir, '-M', 'n', '-Z', tmpdir + '/discrep', '-r', tmpdir
    ])
