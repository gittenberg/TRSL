'''
Created on 31.10.2014

@author: martin
'''
import os.path
from Bio import SeqIO  # @UnresolvedImport

datadir = "C://Users//MJS//Documents//Google Drive//Studium//Master thesis//Data"
sequencefile = "orf_coding.fasta"
sequencefile = os.path.join(datadir, sequencefile)

fasta_sequences = SeqIO.parse(open(sequencefile), 'fasta')
gene_library = {}

for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    gene_library[name] = sequence

for k in sorted(gene_library.keys())[:10]:
    print k
    print gene_library[k]
    print 
'''

# pickle everything
pickle_filename = "gene_library.pkl"
import pickle
import os.path
pickle.dump(gene_library, open(os.path.join(pickle_filename), "wb"))
'''
