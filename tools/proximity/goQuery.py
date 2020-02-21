import synonymParse as sp
import requests, sys
import json

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
        go_data = {}
        for go_term in go_terms:
            go_data[go_term] = {}
            requestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"+str(go_term)

            r = requests.get(requestURL, headers={ "Accept" : "application/json"})
            if not r.ok:
                r.raise_for_status()
                sys.exit()

            term_data = r.json()
            #print(go_term)
            #print('++++++++++++++')
            for go_term_return in term_data['results']: 
                if go_term_return['isObsolete'] == True:
                    print('im here')
                    go_data.pop(go_term)
                else:
                    go_data[go_term]['id'] = go_term_return['id']
                    go_data[go_term]['name'] = go_term_return['name']
                    go_data[go_term]['obsoletion'] = go_term_return['isObsolete']
                    go_data[go_term]['definition'] = go_term_return['definition']['text']
                    go_data[go_term]['usage'] = go_term_return['usage']
                    #print(go_term_return)
                    try:
                        for each_syn in go_term_return['synonyms']:
                            #print(each_syn)
                            if each_syn['type'] == 'exact':
                                go_data[go_term]['synonyms'] = [each_syn['name']]
                            else:
                                continue
                        #print(go_term_return)
                        #print('==========')
                    except KeyError:
                        go_data[go_term]['synonyms'] = ['none to report']
        #{k: v for k, v in go_data.items() if v is not None}
        return go_data
        
    def store_queries(self):
        # What I want is each object that is made from query_gos to be passed to a new dict featuring {initial query : {return results}}
        the_queries = {}
        the_queries[self.term] = self.query_gos()
            
        print(the_queries)
        print('\n')
        print('The size of queries is: '+str(len(the_queries)))
        for k, v in the_queries.items():
            print(k+' : ')
            print('\n')
            for kk, vv in v.items():
                print(kk+' : '+str(vv['name']))
                print(vv['synonyms'])


    def write_list_of_dicts(self):
        pass
'''
############################################################### TEST RANGE / SCRIPT RANGE ##################################################
'''
if __name__ == "__main__":
    my_interest = ['lysis','endolysin','amidase']
    for item in my_interest:
        g = GOQ(term=item)
        #list_of_gos = g.search_gos()
        #print(list_of_gos)
        #g.query_gos()
        g.store_queries()

