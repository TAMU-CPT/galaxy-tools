import synonymParse as sp
import requests, sys

class GOQ:
    """
    GOQ = Gene Ontology Query
    Allows various queries, using QuickGO's REST API.
    --> query : Does a basic query based on input. Retrieves GO terms. Omits obsolete.
    --> expand : For a GO term, fetches data about it.
    """
    def __init__(self, term):
        self.term = term # term that wants to be queried
    
    def search_gos(self):
        ''' searches GO for input term ; queries to a limit of 50 ; returns list of GO terms'''
        requestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/search?query="+str(self.term)+"&limit=50&page=1"

        r = requests.get(requestURL, headers={ "Accept" : "application/json"})

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        q_ret = r.json()
        go_terms = []
        for item in q_ret['results']:
            go_terms.append(item['id'])
        return go_terms
        
    def query_gos(self):
        ''' Searches the return from search gos and finds NON obsolete Go terms; returns dict of GO term | Name | Synonyms '''
        go_terms = self.search_gos()
        for go_term in go_terms:
            requestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"+str(go_term)

            r = requests.get(requestURL, headers={ "Accept" : "application/json"})
            if not r.ok:
                r.raise_for_status()
                sys.exit()

            term_data = r.json()
            print(go_term)
            print('++++++++++++++')
            print(term_data)

        


if __name__ == "__main__":
    g = GOQ(term='lysis')
    #list_of_gos = g.search_gos()
    #print(list_of_gos)
    g.query_gos()

