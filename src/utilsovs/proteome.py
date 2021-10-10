#semantic.py
from Bio import Entrez
from utilsovs.commons import request_API_json

def request_proteomeXchange(data):

    data.url ='http://proteomecentral.proteomexchange.org/cgi/GetDataset?action=search&filterstr='

    data = request_API_json(data)

    return data
