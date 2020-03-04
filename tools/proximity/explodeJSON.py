import json


class explodeJSON:
    def __init__(self, file):
        self.file = file

    def readJSON(self):
        """ returns dictionary object for reading a JSON """
        with open(self.file) as j:
            myObj = json.load(j)

        return myObj

    def explode(self):
        """ Makes a list of each embedded list from the database JSON """

        data = self.readJSON()

        terms = []
        for v in data.values():
            for term in v:
                terms.append(term)

        return terms


if __name__ == "__main__":
    query = []
    filepath = "test-data/"
    filename = "test.json"
    e = explodeJSON(file=filepath + filename)
    data = e.readJSON()
    print(data)
    for k, v in data.items():
        for term in v:
            print(k + ":" + term)  # print global term to synonym / children terms.

    print("++ ========= ++")

    terms = e.explode()
    print(terms)
