#compute.py
from utilsovs.uniprot import request_uniprotKB_API, get_sprot
from utilsovs.globals import PROTEASES, MONOISOTOPIC_AVERAGE_MASS, COLORSCHEME
from utilsovs.commons import init_dict_aa, write_pickle
from utilsovs.prepare import get_aaSprot_freqTable
import math
import copy
import gzip
import re

def compute_combFrags(data):
    data.data = [ data.data[a:b] for a in range(len(data.data)+1) for b in range(len(data.data)+1) if a < b ]
    data.data = [ [ ''.join([ x[0] for x in comb ]), (comb[0][1][0],comb[-1][1][1]) ] for comb in data.data ]
    return data

def compute_splitSeq(data):

    lFrags = re.split(PROTEASES[data.protease][1],data.sequence)

    lIdx = [m.end() for m in re.finditer(PROTEASES[data.protease][1],data.sequence)]

    if len(lFrags) > len(lIdx):
        if len(lFrags) != 1:
            lIdx.append(lIdx[-1])
        else:
            lIdx.append(len(lFrags[0]))

    loFrags = []

    for i, v in enumerate(lFrags):
        if i == 0 and len(v) > 0:
            loFrags.append([v,(1,lIdx[i])])
        elif i == len(lFrags)-1 and len(v) > 0:
            loFrags.append([v,(lIdx[i]+1,len(data.sequence))])
        elif len(v) > 0:
            loFrags.append([v,(lIdx[i-1]+1,lIdx[i])])

    data.data = copy.deepcopy(loFrags)

    return data

def compute_mwFrags(data):
    data.input = copy.deepcopy(data.data)
    data.data = []
    for i in range(len(data.input)):
        wm = MONOISOTOPIC_AVERAGE_MASS['water'][0]
        wa = MONOISOTOPIC_AVERAGE_MASS['water'][0]
        for aa in data.input[i][0]:
            wm = wm + MONOISOTOPIC_AVERAGE_MASS[aa][0]
            wa = wa + MONOISOTOPIC_AVERAGE_MASS[aa][1]
        data.data.append([data.input[i][0],(data.input[i][1][0],data.input[i][1][1]),wm,wa])
    return data

def compute_mw(string):

    wm = MONOISOTOPIC_AVERAGE_MASS['water'][0]
    wa = MONOISOTOPIC_AVERAGE_MASS['water'][0]

    for aa in string:
        wm = wm + MONOISOTOPIC_AVERAGE_MASS[aa][0]
        wa = wa + MONOISOTOPIC_AVERAGE_MASS[aa][1]

    return string, wm, wa

def get_one_match(data):
    data.data = data.id+' '+data.aa+' '+str(data.pos)+' '+data.sequence[data.pos-1]

    if data.aa == data.sequence[data.pos-1]:
        data.data = data.data+' MATCH'
    else:
        data.data = data.data+' MISMATCH'
    return data

def compute_relFreq_seqAlign(data):

    data.data = []

    data.lenSeqs = len(data.input[0])

    for i in range(data.lenSeqs):
        data.data.append({})

        for aa in COLORSCHEME.keys():
            data.data[i][aa] = 0

    for i in range(data.lenSeqs):

        for seq in data.input:
            data.data[i][seq[i]] += 1

    for i in range(data.lenSeqs):

        for aa in data.data[i].keys():
            data.data[i][aa] = data.data[i][aa] / sum(list(data.data[i].values()))

    for i in range(data.lenSeqs):

        for aa in data.data[i].keys():
            try:
                data.data[i][aa] = math.log(data.data[i][aa] / data.aaSprot_freqTable[aa],2)
            except:
                data.data[i][aa] = 0

    return data
