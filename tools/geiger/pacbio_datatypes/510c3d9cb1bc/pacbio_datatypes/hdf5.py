from galaxy.datatypes.binary import H5


class PacBioReads(H5):
    """Class for domain BLAST database files."""
    file_ext = 'pacbioreads'
    allow_datatype_change = False
    composite_type = 'basic'

    def __init__(self, **kwd):
        H5.__init__(self, **kwd)
        self.add_composite_file('parent.bas.h5', is_binary=True)

    def set_peek(self, dataset, is_multi_byte=False):
        """Set the peek and blurb text."""
        if not dataset.dataset.purged:
            dataset.peek = "PacBio Read Set (multiple files)"
            dataset.blurb = "PacBio Read Set (multiple files)"
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'

    def display_peek(self, dataset):
        """Create HTML content, used for displaying peek."""
        try:
            return dataset.peek
        except:
            return "BLAST database (multiple files)"
