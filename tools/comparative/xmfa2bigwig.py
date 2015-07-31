#!/usr/bin/env python
import argparse
import os
import subprocess
import re
import sys
import tempfile
import itertools
from Bio import SeqIO
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def secure_filename(filename):
    r"""Borrowed from :mod:`werkzeug.utils`, under the BSD 3-clause license.
    Pass it a filename and it will return a secure version of it.  This
    filename can then safely be stored on a regular file system and passed
    to :func:`os.path.join`.  The filename returned is an ASCII only string
    for maximum portability.
    On windows systems the function also makes sure that the file is not
    named after one of the special device files.
    >>> secure_filename("My cool movie.mov")
    'My_cool_movie.mov'
    >>> secure_filename("../../../etc/passwd")
    'etc_passwd'
    >>> secure_filename(u'i contain cool \xfcml\xe4uts.txt')
    'i_contain_cool_umlauts.txt'
    The function might return an empty filename.  It's your responsibility
    to ensure that the filename is unique and that you generate random
    filename if the function returned an empty one.
    :param filename: the filename to secure
    """
    PY2 = sys.version_info[0] == 2
    text_type = unicode if PY2 else str

    _filename_ascii_strip_re = re.compile(r'[^A-Za-z0-9_.-]')
    _windows_device_files = ('CON', 'AUX', 'COM1', 'COM2', 'COM3', 'COM4',
                             'LPT1', 'LPT2', 'LPT3', 'PRN', 'NUL')
    if isinstance(filename, text_type):
        from unicodedata import normalize
        filename = normalize('NFKD', filename).encode('ascii', 'ignore')
        if not PY2:
            filename = filename.decode('ascii')
    for sep in os.path.sep, os.path.altsep:
        if sep:
            filename = filename.replace(sep, ' ')
    filename = str(_filename_ascii_strip_re.sub('', '_'.join(
                   filename.split()))).strip('._')

    # on nt a couple of special files are present in each folder.  We
    # have to ensure that the target file is not such a filename.  In
    # this case we prepend an underline
    if os.name == 'nt' and filename and \
       filename.split('.')[0].upper() in _windows_device_files:
        filename = '_' + filename

    return filename


def parse_xmfa(xmfa):
    """Simple XMFA parser until https://github.com/biopython/biopython/pull/544
    """
    current_lcb = []
    current_seq = {}
    for line in xmfa.readlines():
        if line.startswith('#'):
            continue

        if line.strip() == '=':
            if 'id' in current_seq:
                current_lcb.append(current_seq)
                current_seq = {}
            yield current_lcb
            current_lcb = []
        else:
            line = line.strip()
            if line.startswith('>'):
                if 'id' in current_seq:
                    current_lcb.append(current_seq)
                    current_seq = {}
                data = line.strip().split()
                id, loc = data[1].split(':')
                start, end = loc.split('-')
                current_seq = {
                    'rid': '_'.join(data[1:]),
                    'id': id,
                    'start': int(start),
                    'end': int(end),
                    'strand': 1 if data[2] == '+' else -1,
                    'seq': ''
                }
            else:
                current_seq['seq'] += line.strip()


def _percent_identity(a, b):
    """Calculate % identity, ignoring gaps in the host sequence
    """
    match = 0
    mismatch = 0
    for char_a, char_b in zip(list(a), list(b)):
        if char_a == '-':
            continue
        if char_a == char_b:
            match += 1
        else:
            mismatch += 1

    if match + mismatch == 0:
        return 0
    return 100 * float(match) / (match + mismatch)


def _id_tn_dict(sequences):
    """Figure out sequence IDs
    """
    label_convert = {}
    correct_chrom = None
    for i, record in enumerate(SeqIO.parse(sequences, 'fasta')):
        if correct_chrom is None:
            correct_chrom = record.id

        label_convert[str(i + 1)] = {
            'record_id': record.id,
            'len': len(record.seq),
            #'temp': tempfile.NamedTemporaryFile(delete=False)
        }

    return label_convert


def _tempfiles(tn_dict):
    tmpfiles = {}
    for (a, b) in itertools.permutations(tn_dict, 2):
        key = '%s %s' % (a, b)
        tmpfiles[key] = tempfile.NamedTemporaryFile(delete=False, prefix="cpt.xmfa2bw.")
        tmpfiles[key].write("variableStep chrom=%s\n" % tn_dict[a]['record_id'])
    return tmpfiles


def convert_to_bigwig(wig_file, chr_sizes, bw_file):
    #This will be fine under Galaxy, but could use temp folder?
    size_file = "%s-sizes.txt" % (os.path.splitext(bw_file)[0])
    with open(size_file, "w") as out_handle:
        for chrom, size in chr_sizes:
            out_handle.write("%s\t%s\n" % (chrom, size))
    try:
        cl = ["wigToBigWig", wig_file, size_file, bw_file]
        log.debug('CLI: %s', ' '.join(cl))
        subprocess.check_call(cl)
    finally:
        os.remove(wig_file)
        os.remove(size_file)
    return bw_file


def remove_gaps(parent, other):
    # one or more dashes
    m = re.compile('-+')
    fixed_parent = ""
    fixed_other = ""

    last_gap_idx = 0
    # Body
    for gap in m.finditer(parent):
        fixed_parent += parent[last_gap_idx:gap.start()]
        fixed_other += other[last_gap_idx:gap.start()]

        last_gap_idx = gap.end()

    # Tail
    fixed_parent += parent[last_gap_idx:]
    fixed_other += other[last_gap_idx:]

    return fixed_parent, fixed_other


def convert_xmfa_to_gff3(xmfa_file, fasta_genomes, window_size=3):
    label_convert = _id_tn_dict(fasta_genomes)
    tmpfiles = _tempfiles(label_convert)

    try:
        os.makedirs("out")
    except Exception:
        pass

    log.info("Making Comparisons")
    for lcb_idx, lcb in enumerate(parse_xmfa(xmfa_file)):
        log.debug("Processing lcb %s", lcb_idx)
        for (parent, other) in itertools.permutations(lcb, 2):
            local_tmp = tmpfiles['%s %s' % (parent['id'], other['id'])]

            if parent['start'] == 0 and parent['end'] == 0:
                continue

            corrected_parent, corrected_target = remove_gaps(
                parent['seq'],
                other['seq']
            )

            for i in range(1, len(corrected_parent) - 1):
                left_bound = max(0, i - window_size)
                right_bound = i + window_size
                point_pid = _percent_identity(
                    corrected_parent[left_bound:right_bound],
                    corrected_target[left_bound:right_bound]
                )

                local_tmp.write("%s\t%s\n" % (
                    abs(parent['start']) + i,
                    point_pid
                ))

    log.info("Converting files")
    for (pid, oid) in itertools.permutations(label_convert, 2):
        parent = label_convert[pid]
        other = label_convert[oid]

        key = '%s %s' % (pid, oid)
        tmpfiles[key].close()

        bw_file = os.path.join(
            "out", secure_filename(
                '%s - %s.bigwig' % (parent['record_id'], other['record_id']))
        )
        log.debug("Converting %s to %s", key, bw_file)

        # Generate a custom sizes file where all the "TO" genomes are listed as
        # being the same size as the "FROM" (key) genome.
        sizes = [(label_convert[pid]['record_id'], label_convert[pid]['len'])]
        convert_to_bigwig(tmpfiles[key].name, sizes, bw_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert XMFA file to BigWig tracks')
    parser.add_argument('xmfa_file', type=file, help='XMFA file')
    parser.add_argument('fasta_genomes', type=file, help='Fasta genomes')
    parser.add_argument('--version', action='version', version='0.1')
    args = parser.parse_args()

    convert_xmfa_to_gff3(**vars(args))
