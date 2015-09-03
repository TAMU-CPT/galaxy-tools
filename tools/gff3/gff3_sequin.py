import re
from Bio import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


class SequinReportParser(object):

    def gen_qc_feature(start, end, message, strand=0, level='FATAL'):
        return SeqFeature(
            FeatureLocation(start, end, strand=strand),
            qualifiers={
                'note': [message],
                'color': ['#FF0000' if level == 'FATAL' else '#FFFF33'],
            }
        )

    def __init__(self, id='Mordin'):
        self.regex = {
            # Info?
            'DiscRep_ALL': re.compile('^DiscRep_ALL:(?P<ID>[^:]*)::(?P<MESSAGE>.*)'),
            # Info?
            'DiscRep_SUB': re.compile('^DiscRep_SUB:(?P<ID>[^:]*)::(?P<MESSAGE>.*)'),
            # Def. an error
            'DiscRep_Error_Fatal': re.compile('^FATAL: DiscRep_ALL:(?P<ID>[^:]*)::(?P<MESSAGE>.*)'),
        }
        self.id = id

    def parse_discrep(self, discrep):
        gff3_qc_record = SeqRecord(self.id, id=self.id)
        gff3_qc_record.features = []

        with open(discrep, 'r') as handle:
            detailed = False
            # current_id = None
            for line in handle:
                if line.startswith('Detailed Report'):
                    detailed = True

                if detailed:
                    for regex_key in self.regex:
                        m = re.match(self.regex[regex_key], line.strip())
                        if m:
                            if regex_key in ('DiscRep_Erorr_Fatal'):
                                # gen_qc_feature()
                                # # Error
                                # print m.group('ID')
                                pass

if __name__ == '__main__':
    import sys
    srp = SequinReportParser()
    srp.parse_discrep(sys.argv[1])
