import requests
import json
import copy

class GuanineClient(object):

    def __init__(self,
                 url='https://cpt.tamu.edu/guanine/submit',
                 username=None, password=None,
                 api_key='b7c3db072bb74334bd450661c19c6a43',
                 id='9bbd106949b742a5a3961387011e1d5f'):
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
