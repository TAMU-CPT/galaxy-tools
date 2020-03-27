import synonymParse as sp
import requests, sys
import json
import explodeJSON as ej
import argparse


class GOQ:
    """
    GOQ = Gene Ontology Query
    Allows various queries, using QuickGO's REST API.
    --> query : Does a basic query based on input. Retrieves GO terms. Omits obsolete.
    --> expand : For a GO term, fetches data about it.
    """

    def __init__(self, term):
        self.term = term  # term that wants to be queried

    def search_gos(self):
        """ searches GO for input term ; queries to a limit of 50 ; returns list of GO terms"""
        requestURL = (
            "https://www.ebi.ac.uk/QuickGO/services/ontology/go/search?query="
            + str(self.term)
            + "&limit=50&page=1"
        )

        r = requests.get(requestURL, headers={"Accept": "application/json"})

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        q_ret = r.json()
        go_terms = []
        for item in q_ret["results"]:
            go_terms.append(item["id"])
        return go_terms

    def query_gos(self):
        """ Searches the return from search gos and finds NON obsolete Go terms; returns dict of GO term | Name | Synonyms """
        go_terms = self.search_gos()
        go_data = {}
        for go_term in go_terms:
            go_data[go_term] = {}
            requestURL = (
                "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"
                + str(go_term)
            )

            r = requests.get(requestURL, headers={"Accept": "application/json"})
            if not r.ok:
                r.raise_for_status()
                sys.exit()

            term_data = r.json()

            for go_term_return in term_data["results"]:
                if go_term_return["isObsolete"] == True:
                    go_data.pop(go_term)
                else:
                    go_data[go_term]["id"] = go_term_return["id"]
                    go_data[go_term]["name"] = go_term_return["name"]
                    go_data[go_term]["obsoletion"] = go_term_return["isObsolete"]
                    go_data[go_term]["definition"] = go_term_return["definition"][
                        "text"
                    ]
                    go_data[go_term]["usage"] = go_term_return["usage"]

                    try:
                        go_data[go_term]["synonyms"] = []
                        for each_syn in go_term_return["synonyms"]:

                            if each_syn["type"] == "exact":

                                go_data[go_term]["synonyms"] += [each_syn["name"]]
                            else:

                                continue
                        if go_data[go_term]["synonyms"] == []:
                            go_data[go_term]["synonyms"] = [
                                "None"
                            ]  # states that there were no "exact" matching synonyms
                    except KeyError:
                        go_data[go_term]["synonyms"] = [
                            "NoneAtAll"
                        ]  # Means that the return data from GO stated that there were NO synonyms

        return go_data

    def store_queries(self):
        """ Self checking ... can always check this dictionary if the go terms match (value of intial dict key to id value)"""
        the_queries = {}
        the_queries[self.term] = self.query_gos()
        return the_queries

    def write_list_of_dicts(self, filename):

        query_dict = self.store_queries()
        with filename as f:
            # f.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format('query_term','go_id','name','go_definition','usage','synonyms'))
            for search_item, matches in query_dict.items():
                for data in matches.values():
                    f.write(
                        "{}\t{}\t{}\t{}\t{}\t".format(
                            search_item,
                            data["id"],
                            data["name"],
                            data["definition"],
                            data["usage"],
                        )
                    )
                    if len(data["synonyms"]) <= 1:
                        f.write(str(data["synonyms"][0]) + "\n")
                    else:
                        for synonym in data["synonyms"]:
                            if synonym == data["synonyms"][-1]:
                                f.write(
                                    "{}\n".format(synonym)
                                )  # last one doesn't need pipe separator
                            else:
                                f.write("{}|".format(synonym))

                # print('\n')


"""
############################################################### TEST RANGE / SCRIPT RANGE ##################################################
"""
if __name__ == "__main__":

    # Parameters
    parser = argparse.ArgumentParser(
        description="Gene Ontology (GO) term and synonym returns from an input query"
    )
    parser.add_argument(
        "--fileType",
        choices=("database", "new line txt file", "manual insert"),
        default="new line txt file",
        help="help the script decide what kind of file it is going to analyze",
    )  # input file or list
    parser.add_argument(
        "file", type=argparse.FileType("r"), help="Input File"
    )  # Input file --> if testing, use syn.txt
    parser.add_argument(
        "--output",
        dest="output",
        type=argparse.FileType("w"),
        default="go-synonym-results.txt",
        help="Name of the output file",
    )  # output file name

    args = parser.parse_args()

    if args.fileType == "database":
        data = args.file.name
        e = ej.explodeJSON(data)
        terms_to_search = e.explode()
        for term in terms_to_search:
            if term == terms_to_search[0]:
                g = GOQ(term=term)
                g.write_list_of_dicts(filename=args.output)
            else:
                filename = args.output.name
                filename = open(filename, "a+")
                g = GOQ(term=term)
                g.write_list_of_dicts(filename=filename)

    elif args.fileType == "new line txt file":
        terms_to_search = []
        with args.file as f:
            for line in f.readlines():
                line = line.split()[0]
                terms_to_search.append(line)
        for term in terms_to_search:
            if term == terms_to_search[0]:
                g = GOQ(term=term)
                g.write_list_of_dicts(filename=args.output)
            else:
                filename = args.output.name
                filename = open(filename, "a+")
                g = GOQ(term=term)
                g.write_list_of_dicts(filename=filename)

    elif args.fileType == "manual insert":
        pass
