import synonymParse as sp
import requests, sys

class GOQ:
    """
    GOQ = Gene Ontology Query
    Allows various queries, using QuickGO's REST API.
    --> query : Does a basic query based on input. Retrieves GO terms. Omits obsolete.
    --> expand : For a GO term, fetches data about it.
    """
    def __init__(self, file):
        self.file = file
        pass
    
    def query(self):
        ''' searches GO for input term '''
        pass


if __name__ == "__main__":
    query = []
    pass