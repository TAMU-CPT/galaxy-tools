import json

class explodeJSON:

    def __init__(self, file):
        self.file = file

    def readjson(self):
        """ returns dictionary object for reading a JSON """
        with open(self.file) as j:
            myObj = json.load(j) 
        
        return myObj

if __name__ == "__main__":
    query = []
    filepath = 'test-data/'
    filename = 'test.json'
    e = explodeJSON(file=filepath+filename)
    data = e.readjson()
    print(data)
    for k, v in data.items():
        for term in v:
            print(k+ ':' +term) # print global term to synonym / children terms.
