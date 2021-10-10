#semantic.py
from Bio import Entrez
from utilsovs.commons import request_API_json

def request_GlyGen_API(data):

    data.url ='https://api.glygen.org/directsearch/protein/?query='

    data.id = '{"uniprot_canonical_ac":"'+data.id+'"}'

    data = request_API_json(data)

    return data
