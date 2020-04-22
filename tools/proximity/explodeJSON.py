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

### Dictionary Functions
def save_dict_to_json(obj,filename="output.json"):
    with open(filename, "w") as js:
        print("saved {} as json".format(filename))
        json.dump(obj, js, indent=4)


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

    test = {"math":["algebra","calculus"]}
    print(type(test))
    save_dict_to_json(obj=test,filename="test-output.json")