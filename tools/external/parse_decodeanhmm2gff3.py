#!/usr/bin/env python
import sys
import uuid
import argparse
import logging
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC
from bigwig import bigwig_add_header, bigwig_store
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def convert(data=None, bw_i=None, bw_o=None, bw_m=None):
    record = None
    i_data = []
    o_data = []
    m_data = []
    length = None
    count = 1

    if data is None:
        data = sys.stdin

    for idx, line in enumerate(data):
        if line.startswith('>'):
            if record is not None:
                bigwig_store(bw_i, record.id, i_data)
                bigwig_store(bw_o, record.id, o_data)
                bigwig_store(bw_m, record.id, m_data)

            if record is not None and len(record.features) > 0:
                yield record
                # Reset
            i_data = []
            o_data = []
            m_data = []

            header = line.split(' ')[0].strip()[1:]
            record = SeqRecord(
                Seq("ACTG", IUPAC.IUPACUnambiguousDNA),
                id=header
            )

        elif line.startswith('%len'):
            length = int(line[5:].strip())

        elif line.startswith('%pred'):
            pred = line.strip()
            pred = pred[pred.index(': ') + 2:]
            # %pred NB(0): i 1 8, M 9 28, o 29 44
            regions = pred.split(', ')

            tempSub = []

            # Ignore regions that are boring
            if len(regions) == 1:
                continue

            if regions[0].startswith('i'):
                n_dir = 'in'
            else:
                n_dir = 'out'

            if regions[-1].startswith('i'):
                c_dir = 'in'
            else:
                c_dir = 'out'

            # Q7TNJ0  UniProtKB   Chain   1   470 .   .   .   ID=PRO_0000072585;Note=Dendritic cell-specific transmembrane protein

            for region in regions:
                (region_type, start, end) = region.strip().split(' ')
                qualifiers = {
                    'source': 'TMHMM',
                    # 'evidence': 'ECO:0000255',
                }
                # Q7TNJ0  UniProtKB   Topological domain  1   33  .   .   .   Note=Cytoplasmic;evidence=ECO:0000255
                # Q7TNJ0  UniProtKB   Transmembrane   34  54  .   .   .   Note=Helical;evidence=ECO:0000255
                # Q7TNJ0  UniProtKB   Topological domain  55  55  .   .   .   Note=Extracellular;evidence=ECO:0000255
                if region_type == 'i':
                    qualifiers.update({
                        'Note': 'Cytoplasmic',
                    })
                elif region_type == 'M':
                    qualifiers.update({
                        'Note': 'Helical',
                    })
                elif region_type == 'o':
                    qualifiers.update({
                        'Note': 'Extracellular',
                    })

                sub_feat = SeqFeature(
                    FeatureLocation(int(start) - 1, int(end)),
                    type="Transmembrane" if region_type == 'M' else "Topological domain",
                    strand=1,
                    qualifiers=qualifiers,
                )
                tempSub.append(sub_feat)

            feature = SeqFeature(
                FeatureLocation(1, length),
                type="Chain",
                strand=1,
                qualifiers={
                    'ID': 'tmhmm_tmd_%s-%s' % (count, str(uuid.uuid4())),
                    'Description': 'Transmembrane protein',
                    'Note': 'Transmembrane protein - N %s C %s' % (n_dir, c_dir),
                    'Target': header,
                },
                sub_features = tempSub
            )
            count += 1
            

            record.features.append(feature)
        else:
            if record:
                parts = line.split()
                if len(parts) != 5:
                    continue

                if parts[0] == '#':
                    continue

                # ['#', 'i', 'O', 'o', 'M']
                # in, out, ?, membrane potential
                i_data.append(float(parts[1]))
                o_data.append(float(parts[3]))
                m_data.append(float(parts[4]))

    bigwig_store(bw_i, record.id, i_data)
    bigwig_store(bw_o, record.id, o_data)
    bigwig_store(bw_m, record.id, m_data)
    yield record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process TMHMM outputs in GFF3, BigWig')
    parser.add_argument('--bw_i', default='tmhmm_i.wig')
    parser.add_argument('--bw_o', default='tmhmm_o.wig')
    parser.add_argument('--bw_m', default='tmhmm_m.wig')
    args = parser.parse_args()

    bw_i = open(args.bw_i, 'w')
    bw_o = open(args.bw_o, 'w')
    bw_m = open(args.bw_m, 'w')

    bigwig_add_header(bw_i, 'i', name='TMHMM')
    bigwig_add_header(bw_o, 'o', name='TMHMM')
    bigwig_add_header(bw_m, 'm', name='TMHMM')

    for sequence in convert(None, bw_i, bw_o, bw_m):
        GFF.write([sequence], sys.stdout)

    bw_i.close()
    bw_o.close()
    bw_m.close()
