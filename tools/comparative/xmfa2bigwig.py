#!/usr/bin/env python
import argparse
import os
import subprocess
import re
import sys
import tempfile
from Bio import SeqIO
from xmfa import parse_xmfa, percent_identity


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


def _id_tn_dict(sequences):
    """Figure out sequence IDs
    """
    label_convert = {}
    correct_chrom = None
    for i, record in enumerate(SeqIO.parse(sequences, 'fasta')):
        print i, record
        if correct_chrom is None:
            correct_chrom = record.id

        label_convert[str(i + 1)] = {'record_id': record.id,
                                     'len': len(record.seq),
                                     'temp': tempfile.NamedTemporaryFile(delete=False)}

        label_convert[str(i + 1)]['temp'].write("variableStep chrom=%s\n" % (correct_chrom or record.id, ))
    return label_convert


def convert_to_bigwig(wig_file, chr_sizes, bw_file):
    # This will be fine under Galaxy, but could use temp folder?
    size_file = "%s-sizes.txt" % (os.path.splitext(bw_file)[0])
    with open(size_file, "w") as out_handle:
        for chrom, size in chr_sizes:
            out_handle.write("%s\t%s\n" % (chrom, size))
    try:
        cl = ["wigToBigWig", wig_file, size_file, bw_file]
        print ' '.join(cl)
        subprocess.check_call(cl)
    finally:
        pass
        os.remove(wig_file)
        os.remove(size_file)
    return bw_file


def remove_gaps(parent, others):
    # one or more dashes
    m = re.compile('-+')
    fixed_parent = ""
    fixed_others = ["" for x in others]

    last_gap_idx = 0
    # Body
    for gap in m.finditer(parent):
        fixed_parent += parent[last_gap_idx:gap.start()]
        for i, other in enumerate(others):
            fixed_others[i] += other[last_gap_idx:gap.start()]
        last_gap_idx = gap.end()

    # Tail
    fixed_parent += parent[last_gap_idx:]
    for i, other in enumerate(others):
        fixed_others[i] += other[last_gap_idx:]

    return fixed_parent, fixed_others


def convert_xmfa_to_gff3(xmfa_file, fasta_genomes, window_size=3, relative_to='1'):
    label_convert = _id_tn_dict(fasta_genomes)
    try:
        os.makedirs("out")
    except Exception:
        pass

    for lcb_idx, lcb in enumerate(parse_xmfa(xmfa_file)):
        ids = [seq['id'] for seq in lcb]

        # Doesn't match part of our sequence
        if relative_to not in ids:
            continue

        # Skip sequences that are JUST our "relative_to" genome
        if len(ids) == 1:
            continue

        parent = [seq for seq in lcb if seq['id'] == relative_to][0]
        others = [seq for seq in lcb if seq['id'] != relative_to]

        if parent['start'] == 0 and parent['end'] == 0:
            continue

        corrected_parent, corrected_targets = remove_gaps(parent['seq'],
                                                          [other['seq'] for other in others])
        # Update the parent/others with corrected sequences
        parent['corrected'] = corrected_parent
        for i, target in enumerate(corrected_targets):
            others[i]['corrected'] = target

        for i in range(1, len(corrected_parent) - 1):
            for other in others:
                left_bound = max(0, i - window_size)
                right_bound = i + window_size
                point_pid = percent_identity(
                    parent['corrected'][left_bound:right_bound],
                    other['corrected'][left_bound:right_bound]
                )

                label_convert[other['id']]['temp'].write("%s\t%s\n" % (
                    abs(parent['start']) + i,
                    point_pid
                ))

    for key in label_convert.keys():
        # Ignore self-self
        if key == relative_to:
            continue

        other = label_convert[key]
        other['temp'].close()
        sizes = [(label_convert[relative_to]['record_id'], label_convert[relative_to]['len'])]
        bw_file = os.path.join("out", secure_filename(other['record_id'] + '.bigwig'))

        convert_to_bigwig(label_convert[key]['temp'].name, sizes, bw_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert XMFA file to BigWig tracks')
    parser.add_argument('xmfa_file', type=argparse.FileType("r"), help='XMFA file')
    parser.add_argument('fasta_genomes', type=argparse.FileType("r"), help='Fasta genomes')
    parser.add_argument('--window_size', type=int, help='Window Size', default=3)
    parser.add_argument('--version', action='version', version='0.1')
    args = parser.parse_args()

    convert_xmfa_to_gff3(**vars(args))
