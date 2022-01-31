# Utilsovs - 0.9.3

Utils derived from the [O-GlcNAc Database](https://www.oglcnac.mcw.edu/) code source.

Please report any bugs or incompatibilities.

If you use *utilsovs* in your academic work, please cite:

Malard F, Wulff-Fuentes E, Berendt R, Didier G and Olivier-Van Stichelen S. **Automatization and self-maintenance of the O-GlcNAcome catalogue:
A Smart Scientific Database**. *Database*, Volume 2021, (2021).

## Install

```python
pip3 install utilsovs-pkg
```

Test install with ```pytest``` from the package root directory.

## Content

The package utilsovs contains:

- API wrappers - Proteins from UniProtKB ID ([UniProtKB](https://www.uniprot.org/), [GlyGen](https://www.glygen.org/), [The *O*-GlcNAc Database](https://www.oglcnac.mcw.edu/))
- API wrappers - Literature from PMID ([MedLine/PubMed](https://pubmed.ncbi.nlm.nih.gov/), [Semantic Scholar](https://www.semanticscholar.org/), [ProteomeXchange](http://www.proteomexchange.org/))
- Protein digestion tool: full and partial digestion and MW calculation (monoisotopic, average mass)
- Calculation of log2(odds) from alignment file and generation of sequence logo
- Match residuePosition on sequence fetched from UniProtKB to validate datasets
- Convert PDF to Text using wrappers and repair/clean
- Miscellaneous functions

### API wrappers - Proteins from UniProtKB ID

```python
from utilsovs import *

# Fetch UniProtKB Proteins REST API (@data.url)
data = fetch_one_UniProtKB('P08047',filepath='out.json',pprint=False)

# Fetch The O-GlcNAc Database Proteins REST API (@data.url)
data = fetch_one_oglcnacDB('P08047',filepath='out.json',pprint=False)

# Fetch RESTful Glygen webservice-based APIs (@data.url)
data = fetch_one_GlyGen('P08047',filepath='out.json',pprint=False)

# data is an class instance. To print the data of interest:
print (data.data)

```

### API wrappers - Literature from PubMed IDentifier (PMID)

```python
from utilsovs import *

# Fetch MedLine/PubMed API using Entrez.efetch (@data.url)
data = fetch_one_PubMed('33479245',db="pubmed",filepath='out.json',pprint=False)

# Fetch Semantic Scholar API (@data.url)
data = fetch_one_SemanticScholar('33479245',filepath='out.json',pprint=False)

# Fetch proteomeXchange using GET search request (@data.url)
data = fetch_one_proteomeXchange('29351928',filepath='out.json',pprint=False)

# data is an class instance. To print the data of interest:
print (data.data)

```

### Compute - Digest protein, match residuePosition on sequence or calculate log2(odds) from alignment file and draw consensus sequence logo

```python
from utilsovs import *

# Full digestion of a UniProtKB ID protein sequence: [ ['PEPTIDE',(start,end),mw_monoisotopic,mw_average], ... ]
data = compute_one_fullDigest('P13693','Trypsin',filepath='out.json')

# Partial digestion of a UniProtKB ID protein sequence: [ ['PEPTIDE',(start,end),mw_monoisotopic,mw_average], ... ]
# All possible combinations of adjacent fragments are generated
data = compute_one_partialDigest('P13693','Trypsin',filepath='out.json')

# Match residuePosition with UniProtKB ID protein sequence
data = compute_match_aaSeq('P13693','D6',filepath='out.json')

# Compute log2odds from alignment file - Input for draw_one_seqLogo()
data = compute_aln_log2odds('align.aln',organism='HUMAN',filepath='out.json')

# Draw sequence logo from compute_aln_log2odds output file
# See https://logomaker.readthedocs.io/en/latest/implementation.html
# Edit logomaker config in src/ultilsovs_draw.py
draw_one_seqLogo('compute_aln_log2odds.json',filepath='out.png',showplot=False,center_values=False)

# data is an class instance. To print the data of interest:
print (data.data)

```

### Text Processing

```python
from utilsovs import *

# PDF to Text conversion using GNU pdftotext (Linux/Mac) or Tika (Windows) and text repair + cleaning.
data = pdf_one_pdf2text('test.pdf',filepath='out.dat',clean=True)

# data is an class instance. To print the data of interest:
print (data.data)

```

### Miscellaneous standalone functions

Functions below return Python objects or variables.

```python
from utilsovs import *

# Show list of proteases for digest utils
show_proteases()

# Return protein sequence from UniProtKB ID
get_one_sequence('P13693',filepath='out.dat')

# Compute MW of a peptide and return [string,mw_monoisotopic,mw_average]
compute_one_MW('EWENMR',filepath='out.json')

#Compute amino-acids frequency table for a given organism from uniprot_sprot.fasta.gz
get_one_freqAAdict(organism='HUMAN',filepath='out.json')

#Clear all data in utilsovs cache
clearCache()


```
