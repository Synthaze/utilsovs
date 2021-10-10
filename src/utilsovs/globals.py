#globals.py

URL_SPROT = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz'

PROTEASES = {
### This is unstested for N-Ter cut
'Trypsin':['bovine',r'(?<=[KR])(?<=[A-Z])','Cut C-Ter K, R'],
'Chymotrypsin':['bovine',r'(?<=[YFW])(?<=[A-Z])','Cut C-Ter Y, F, W'],
'Arg-C':['mouse-submaxillary-gland',r'(?<=[R])(?<=[A-Z])','Cut C-Ter R'],
'Glu-C':['Staphylococcus-aureus',r'(?<=[E])(?<=[A-Z])','Cut C-Ter E'],
'Lys-C':['Lysobacter-enzymogenes',r'(?<=[K])(?<=[A-Z])','Cut C-Ter K'],
'Pepsin':['porcine',r'(?<=[GAVLIPMYFW])(?<=[A-Z])','Cut C-Ter G, A, V, L, I, P, M, Y, F, W'],
'Thermolysin':['Bacillus-thermo-proteolyticus',r'(?<=[LFIVMA])(?<=[A-Z])','Cut C-Ter L, F, I, V, M, A'],
'Elastase':['porcine',r'(?<=[AVSGLI])(?<=[A-Z])','Cut C-Ter A, V, S, G, L, I'],
'DUMMY':['DUMMY','(?<=[C])(?<=[A-Z])','For testing purpose only'],
#'ADD_ONE':['SPECIES','REGEX','DESCRIPTION'],
}

MONOISOTOPIC_AVERAGE_MASS = {
        'A':[71.03711,71.0788],
        'R':[156.10111,156.1875],
        'N':[114.04293,114.1038],
        'D':[115.02694,115.0886],
        'C':[103.00919,103.1388],
        'E':[129.04259,129.1155],
        'Q':[128.05858,128.1307],
        'G':[57.02146,57.0519],
        'H':[137.05891,137.1411],
        'I':[113.08406,113.1594],
        'L':[113.08406,113.1594],
        'K':[128.09496,128.1741],
        'M':[131.04049,131.1926],
        'F':[147.06841,147.1766],
        'P':[97.05276,97.1167],
        'S':[87.03203,87.0782],
        'T':[101.04768,101.1051],
        'W':[186.07931,186.2132],
        'Y':[163.06333,163.1760],
        'V':[99.06841,99.1326],
        'U':[150.953636,150.0388],
        'O':[237.147727,237.3018],
        'water':[18.01056,18.01524],
        'oglcnac':[203.0794,203.1950],
        'phospho':[79.9663,79.9799]
}

COLORSCHEME = {
    'R': 'blue',
    'K': 'blue',
    'D': 'blue',
    'E': 'blue',
    'N': 'blue',
    'Q': 'blue',
    'S': 'darkgreen',
    'G': 'darkgreen',
    'H': 'darkgreen',
    'T': 'darkgreen',
    'A': 'darkgreen',
    'P': 'darkgreen',
    'Y': 'black',
    'V': 'black',
    'M': 'black',
    'C': 'black',
    'L': 'black',
    'F': 'black',
    'I': 'black',
    'W': 'black',
}
