import requests
import json
import copy

class GuanineClient(object):

    def __init__(self, url='http://localhost:8000/gses/submit',
                 username=None, password=None, api_key=None, id=None):
        self.url = url
        self.default_data = {
            'id': id,
        }
        if api_key is not None:
            self.default_data.update({
                'api_key': api_key
            })
        else:
            self.default_data.update({
                'username': username,
                'password': password,
            })

    def _submit(self, results):
        ldata = copy.copy(self.default_data)
        ldata.update({
            'data': results
        })
        r = requests.post(
            self.url,
            data=json.dumps(ldata)
        )
        return r.text

    def submit(self, user, assessment, score):
        return self._submit([{
            'user': user,
            'assessment': assessment,
            'score': score
        }])
