#semantic.py
from Bio import Entrez
from utilsovs.commons import request_API_json


def request_oglcnacDB_API(data):

    data.url ='https://www.oglcnac.mcw.edu/api/v1/proteins/?query_protein='

    data = request_API_json(data)

    return data
