#uniprot.py
from utilsovs.commons import request_API_json
from utilsovs.globals import URL_SPROT
import pathlib
import wget
import os

def request_uniprotKB_API(data):

    data.url = "https://www.ebi.ac.uk/proteins/api/proteins/"

    data = request_API_json(data)

    return data

def get_sprot(uniprot_sprot_path,download=False):

    if download == True:
        print ('Cache recently cleared or fresh install: downloading uniprot_sprot.fasta.gz from ftp://ftp.uniprot.org')
        wget.download(URL_SPROT, out=uniprot_sprot_path)
    else:
        print ('Reading %s' % uniprot_sprot_path)

def get_one_sequence_auto(data):

    id_fasta_path = str(pathlib.Path(__file__).parent.absolute())+'/fasta/'+data.id+'.fasta'

    #if not os.path.exists(id_fasta_path):
    data.sequence = request_uniprotKB_API(data).data['sequence']['sequence']
        #with open(id_fasta_path,'w') as msg:
        #    msg.write('>'+data.id+'\n')
        #    msg.write(data.sequence)
    #else:
    #    print ('Read cached sequence for %s from %s' % (data.id,id_fasta_path))
    #    with open(id_fasta_path,'r') as msg:
    #       data.sequence = msg.read().splitlines()[1]

    return data
