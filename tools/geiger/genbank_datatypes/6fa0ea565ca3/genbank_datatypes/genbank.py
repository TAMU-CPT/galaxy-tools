"""
CPT Genbank Class
"""
import logging
import re
from galaxy.datatypes.data import get_file_peek
from galaxy.datatypes.data import Text

log = logging.getLogger(__name__)


class GenBank(Text):
    """ Sequence Read Archive (SRA) """
    file_ext = 'gbk'

    def __init__(self, **kwd):
        Text.__init__(self, **kwd)

    def sniff(self, filename):
        """
        The first 12 characters ought to be 'LOCUS       '
        """
        try:
            header = open(filename).read(8)
            if header == 'LOCUS       ':
                return True
            else:
                return False
        except:
            return False

    def set_peek(self, dataset, is_multi_byte=False):
        """Set the peek/blurb text"""
        if not dataset.dataset.purged:
            dataset.peek = get_file_peek(dataset.file_name,
                                         is_multi_byte=is_multi_byte)

            # Blurb creation
            p = re.compile('LOCUS       ([^ ]*)')
            fh = file(dataset.file_name)
            locus_tags = []
            for line in fh:
                if p.match(line):
                    m = p.match(line)
                    locus_tags.append(m.group(1))
            dataset.blurb = "Genbank file with %s genomes: %s", (
                len(locus_tags), ', '.join(locus_tags))
        else:
            dataset.peek = "file does not exist"
            dataset.blurb = "file purged from disk"
