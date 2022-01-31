from utilsovs.models import proteinData, pmidData, digestData, genData, matchData, alignData, dataData
from utilsovs.commons import write_file, pretty_json, read_file, write_json, read_json
from utilsovs.uniprot import request_uniprotKB_API, get_one_sequence_auto
from utilsovs.oglcnacDB import request_oglcnacDB_API
from utilsovs.glygen import request_GlyGen_API
from utilsovs.pubmed import request_MedLine_API, one_pdf2text, clean_pdf2text, search_MedLine_API, download_PubMed_OA
from utilsovs.semantic import request_SemanticScholar_API
from utilsovs.proteome import request_proteomeXchange
from utilsovs.prepare import get_aaSprot_freqTable
from utilsovs.compute import compute_splitSeq, compute_combFrags, compute_mwFrags, get_one_match, compute_relFreq_seqAlign, compute_mw
from utilsovs.draw import draw_seqLogo
from utilsovs.globals import PROTEASES
import pandas as pd
import pathlib
import inspect
import glob
import sys
import os

#### Fetch protein & PMID info from major Database

## Proteins; id must be UniProtKB identifier

# Fetch UniProtKB Proteins REST API (@data.url)
def fetch_one_UniProtKB(id,filepath=None,pprint=True):
    data = proteinData(id,filepath,inspect.stack()[0][3])
    data = request_uniprotKB_API(data)
    if pprint == True:
        pretty_json(data)
    write_json(data)
    return data

# Fetch The O-GlcNAc Database Proteins REST API (@data.url)
def fetch_one_oglcnacDB(id,filepath=None,pprint=True):
    data = proteinData(id,filepath,inspect.stack()[0][3])
    data = request_oglcnacDB_API(data)
    if pprint == True:
        pretty_json(data)
    write_json(data)
    return data

# Fetch RESTful Glygen webservice-based APIs (@data.url)
def fetch_one_GlyGen(id,filepath=None,pprint=True):
    data = proteinData(id,filepath,inspect.stack()[0][3])
    data = request_GlyGen_API(data)
    if pprint == True:
        pretty_json(data)
    write_json(data)
    return data

## PMIDs; id must be PubMed IDentifier (PMID)

# Fetch MedLine/PubMed API using Entrez.efetch (@data.url)
def fetch_one_PubMed(id,db="pubmed",filepath=None,pprint=True):
    data = pmidData(id,filepath,inspect.stack()[0][3])
    data = request_MedLine_API(data,db)
    if pprint == True:
        pretty_json(data)
    write_json(data)
    return data

# Fetch Semantic Scholar API (@data.url)
def fetch_one_SemanticScholar(id,filepath=None,pprint=True):
    data = pmidData(id,filepath,inspect.stack()[0][3])
    data = request_SemanticScholar_API(data)
    if pprint == True:
        pretty_json(data)
    write_json(data)
    return data

# Fetch proteomeXchange using GET search request (@data.url)
def fetch_one_proteomeXchange(id,filepath=None,pprint=True):
    data = pmidData(id,filepath,inspect.stack()[0][3])
    data = request_proteomeXchange(data)
    if pprint == True:
        pretty_json(data)
    write_json(data)
    return data


# Search MedLine/PubMed API using Entrez.efetch (@data.url)
def search_one_PubMed(id, db="pubmed", filepath=None, pprint=True):
    data = pmidData(id,filepath,inspect.stack()[0][3])
    data = search_MedLine_API(data,db)
    if pprint == True:
        pretty_json(data)
    write_json(data)
    return data


# Download from pubmed open access (@data.url)
def download_one_PubMed(id, filepath=None):
    data = pmidData(id,filepath,inspect.stack()[0][3])
    data = download_PubMed_OA(data)
    return data

### Compute utils

## Proteins

# Full digestion of a UniProtKB ID protein sequence: [ ['PEPTIDE',(start,end),mw_monoisotopic,mw_average], ... ]
def compute_one_fullDigest(id,protease,filepath=None):
    data = digestData(id,protease,filepath,inspect.stack()[0][3])
    data = get_one_sequence_auto(data)
    data = compute_splitSeq(data)
    data = compute_mwFrags(data)
    write_json(data)
    return data

# Partial digestion of a UniProtKB ID protein sequence: [ ['PEPTIDE',(start,end),mw_monoisotopic,mw_average], ... ] - All possible combinations are generated
def compute_one_partialDigest(id,protease,filepath=None):
    data = digestData(id,protease,filepath,inspect.stack()[0][3])
    data = get_one_sequence_auto(data)
    data = compute_splitSeq(data)
    data = compute_combFrags(data)
    data = compute_mwFrags(data)
    write_json(data)
    return data

# Match residuePosition with UniProtKB ID protein sequence
def compute_match_aaSeq(id,res,filepath=None):
    data = matchData(id,filepath,res,inspect.stack()[0][3])
    data = get_one_sequence_auto(data)
    data = get_one_match(data)
    write_file(data)
    return data

# Compute log2odds for alignment file - Input for draw_one_seqLogo()
def compute_aln_log2odds(alnpath,organism='HUMAN',filepath=None):
    data = alignData(filepath,alnpath,inspect.stack()[0][3])
    data.organism = organism
    data = get_aaSprot_freqTable(data)
    data = compute_relFreq_seqAlign(data)
    write_json(data)
    return data

## PMIDs

### Draw utils

# Draw sequence logo from compute_aln_log2odds output file - See https://logomaker.readthedocs.io/en/latest/implementation.html for more options and edit the corresponding function in src/ultilsovs_draw.py
def draw_one_seqLogo(infile,filepath=None,showplot=True,center_values=False):
    data = dataData(filepath,inspect.stack()[0][3])
    data.df = pd.read_json(infile)
    draw_seqLogo(data,showplot,center_values)
    return data

# PDF to Text conversion using GNU pdftotext (Linux/Mac) or Tika (Windows) and text repair + cleaning.
def pdf_one_pdf2text(pdfpath,filepath=None,clean=False):
    data = pmidData('DUMMY',filepath,inspect.stack()[0][3])
    data.input = pdfpath
    data = one_pdf2text(data)
    if clean == True:
        data = clean_pdf2text(data)
    write_file(data)
    return data


### Miscellaneous standalone functions

# Show list of proteases for digest utils
def show_proteases():
    print("{:<15} {:<30} {:<30} {:<15}".format('PROTEASE','SPECIES','REGEX','DESCRIPTION'))
    for protease, regex in PROTEASES.items():
        print("{:<15} {:<30} {:<30} {:<15}".format(protease,regex[0],regex[1],regex[2]))
    print (' ')
    print ('The PROTEASES dictionary can be edited at src/utilsovs/globals.py')

# Return protein sequence from UniProtKB ID
def get_one_sequence(id,filepath=None):
    data = proteinData(id,filepath,inspect.stack()[0][3])
    data.data = get_one_sequence_auto(data).sequence
    write_json(data)
    return data.data

# Compute MW of a peptide and return [string,mw_monoisotopic,mw_average]
def compute_one_MW(string,filepath=None):
    data = genData(filepath,[],inspect.stack()[0][3])
    data.data = compute_mw(string)
    write_json(data)
    return data.data

# Compute amino-acids frequency table for a given organism from uniprot_sprot.fasta.gz
def get_one_freqAAdict(organism='HUMAN',filepath=None):
    data = dataData(filepath,inspect.stack()[0][3])
    data.organism = organism
    data = get_aaSprot_freqTable(data)
    data.data = data.aaSprot_freqTable
    write_json(data)
    return data

# Clear all data in utilsovs cache
def clearCache():
    for fname in glob.glob(str(pathlib.Path(__file__).parent.absolute())+'/fasta/*fasta'):
        os.remove(fname)
        print ('fasta/%s deleted...' % fname)
    for fname in glob.glob(str(pathlib.Path(__file__).parent.absolute())+'/fasta/*gz'):
        os.remove(fname)
        print ('fasta/%s deleted...' % fname)
    for fname in glob.glob(str(pathlib.Path(__file__).parent.absolute())+'/tables/*pickle'):
        os.remove(fname)
        print ('tables/%s deleted...' % fname)
    for fname in glob.glob(str(pathlib.Path(__file__).parent.absolute())+'/../../tests/*json'):
        os.remove(fname)
        print ('tests/%s deleted...' % fname)
    for fname in glob.glob(str(pathlib.Path(__file__).parent.absolute())+'/../../tests/*dat'):
        os.remove(fname)
        print ('tests/%s deleted...' % fname)
    for fname in glob.glob(str(pathlib.Path(__file__).parent.absolute())+'/../../tests/*png'):
        os.remove(fname)
        print ('tests/%s deleted...' % fname)
