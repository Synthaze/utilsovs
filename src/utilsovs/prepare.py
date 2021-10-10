#prepare.py
from utilsovs.uniprot import get_sprot
from utilsovs.commons import init_dict_aa, write_pickle, read_pickle
from utilsovs.globals import COLORSCHEME
import pathlib
import gzip
import os

def get_aaSprot_freqTable(data,erase=False):

    data.path_aaSprot_freqTable = str(pathlib.Path(__file__).parent.absolute())+'/tables/aaSprot_freqTable_'+data.organism+'.pickle'

    if not os.path.exists(data.path_aaSprot_freqTable) or erase == True:
        data = prepare_aaSprot_occurences(data)
        prepare_aaSprot_frequencies(data)

    data.aaSprot_freqTable = read_pickle(data.path_aaSprot_freqTable)

    return data

def prepare_aaSprot_occurences(data):

    data.path_sprot = str(pathlib.Path(__file__).parent.absolute())+'/fasta/uniprot_sprot.fasta.gz'

    if not os.path.exists(data.path_sprot):
        download = True
    else:
        download = False

    get_sprot(data.path_sprot,download=download)

    data.occurences = init_dict_aa()

    with gzip.open(data.path_sprot, 'rt') as f:
        for l in f:
            ft_l = l.strip()
            if 'sp' in ft_l:
                test_organism = ft_l.split("|")[2].split()[0].split("_")[1]

                if test_organism == data.organism:
                    read = True
                else:
                    read = False

            if read == True and 'sp' not in ft_l:
                for aa in ft_l:
                    if aa in COLORSCHEME.keys():
                        data.occurences[aa] += 1
        f.close()

    return data

def prepare_aaSprot_frequencies(data):
    data.frequencies = {}
    data.total_aa = sum(list(data.occurences.values()))

    for aa in data.occurences.keys():
        data.frequencies[aa] = data.occurences[aa] / data.total_aa

    data.path_aaSprot_freqTable = str(pathlib.Path(__file__).parent.absolute())+'/tables/aaSprot_freqTable_'+data.organism+'.pickle'

    write_pickle(data.path_aaSprot_freqTable,data.frequencies)
