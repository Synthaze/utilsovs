#pubmed.py
from tika import parser
from Bio import Entrez
from utilsovs.commons import read_file
import platform
import os
import re
import json
import requests


def request_MedLine_API(data, db):

    Entrez.email = 'email@server.com'

    data.data = Entrez.efetch(db=db, id=data.id, retmode="xml")

    data.data = Entrez.read(data.data)

    return data


def search_MedLine_API(data, db):

    Entrez.email = 'email@server.com'

    data.data = Entrez.esearch(db=db,
                               sort='relevance',
                               retmax='100000',
                               retmode='json',
                               term=data.id)

    data.data = json.load(data.data)['esearchresult']

    return data


def altsearch_MedLine_API(request, db):

    Entrez.email = 'email@server.com'

    data = Entrez.esearch(db=db,
                          sort='relevance',
                          retmax='100000',
                          retmode='json',
                          term=request)

    data = json.load(data)['esearchresult']

    return data


def download_PubMed_OA(data):

    data.url = 'https://pubmed.ncbi.nlm.nih.gov/' + data.id

    print(data.url)

    r = requests.get(data.url).text

    try:
        data.url = 'https://www.ncbi.nlm.nih.gov/pmc/articles/'

        for x in r.split('/pmc/articles/'):
            print(x.split('"')[0])

        data.url = data.url + x.split('"')[0]

        print(data.url)

    except:
        print(False)

    return data


def one_pdf2text(data):
    if platform.system() == 'Windows':
        print ('/!\\ Windows detected: using tika')
        data.data = parser.from_file(data.input)['content']
    else:
        print ('/!\\ Linux/Mac detected: using pdftotext')
        os.system('pdftotext %s tmp.dat' % (data.input))
        data.data = read_file('tmp.dat')
    return data

def clean_pdf2text(data):

    data.data = re.sub(r'\n',' ',data.data)

    data.data = re.sub(r'\s+',' ',data.data)

    data.data = data.data.split(".")

    for i, line in enumerate(data.data):

        data.data[i] = re.sub(r'(\[[^\[\]]+\]|\([^\(\)]+\))','',line).strip()

    data.data = [x.strip() for x in data.data if x.strip()]

    data.data = " . ".join(data.data)

    data.data = re.sub(r'[^A-Za-z0-9\s.]',' ',data.data)

    data.data = re.sub(r'\s+',' ',data.data)

    return data
