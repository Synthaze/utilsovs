#commons.py
from utilsovs.globals import COLORSCHEME
from pygments import highlight, lexers, formatters
import requests
import pickle
import time
import json
import sys
import re

def request_API_json(data):
    data.id = data.id.split('-')[0] if data.id.endswith('-0') else data.id
    data.url = data.url + data.id

    n = 0

    while n < 5:
        try:
            r = requests.get(data.url, headers={"Accept": "application/json"},verify=False)
            break
        except:
            time.sleep(5)
            n += 1

    try:
        data.data = json.loads(r.text)
    except:
        data.data = ('Error - Can not be found at %s' % data.url)
    return data

def pretty_json(data):
    print (highlight(json.dumps(data.data,sort_keys=True,indent=4), lexers.JsonLexer(), formatters.TerminalFormatter()))
    return data

def write_json(data):
    if data.filepath != None:
        try:
            with open(data.filepath,'w') as msg:
                json.dump(data.data,msg)
            print ('JSON File %s successfully written to disk by %s routine' % (data.filepath,data.func))
        except:
            print ('Program failed, your request may not be processed. Terminating' % (data.filepath,data.func))
    else:
        pass

def read_json(f):
    with open(f) as msg:
        data = json.load(msg)
    return data

def write_file(data):
    if data.filepath != None:
        with open(data.filepath,'w') as msg:
            msg.write(data.data)
        print ('File %s successfully written to disk by %s routine' % (data.filepath,data.func))
    else:
        pass

def read_file(f):
    with open(f, 'r') as fp:
        data = fp.read()
    return data

def init_dict_aa():
    d = {}
    for aa in COLORSCHEME.keys():
        d[aa] = 0
    return d

def write_pickle(f,c):
    with open(f, 'wb') as msg:
        pickle.dump(c,msg)

def read_pickle(f):
    with open(f, 'rb') as msg:
        c = pickle.load(msg)
    return c
