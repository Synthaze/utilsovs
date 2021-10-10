#commons.py
from utilsovs.globals import COLORSCHEME
from pygments import highlight, lexers, formatters
import requests
import pickle
import json
import sys
import re

def request_API_json(data):
    data.url = data.url + data.id
    r = requests.get(data.url, headers={"Accept": "application/json"},verify=False)
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

def control_PMID(id):
    if re.match(r"^[0-9]+$",id):
        print ("PMID %s is valid" % id)
    else:
        print ("PMID %s is not a valid PMID, sys.exit() now..." % id)
        sys.exit()

def control_UniProtKB(id):
    if re.match(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$",id):
        print ("UniProtKB ID %s is valid" % id)
    else:
        print ("UniProtKB ID %s is not a valid UniProtKB ID, sys.exit() now..." % id)
        sys.exit()

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
