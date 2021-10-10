#semantic.py
from utilsovs.commons import request_API_json

def request_SemanticScholar_API(data):

    data.url ='https://api.semanticscholar.org/v1/paper/PMID:'

    data = request_API_json(data)

    return data
