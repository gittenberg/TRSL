'''
Created on 07.01.2015

@author: MJS

TODO: convert to ipython notebook!!
'''
import os.path
import csv
import cPickle as pkl

from Bio import SeqIO

datadir = "data"
# Transcriptome: http://www.ncbi.nlm.nih.gov/pubmed/19581875
transcriptomefile = r"nbt.1551-S2.csv"
# Exome: http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/
sequencefile = "orf_coding.fasta"
transcriptomefilepath = os.path.join(datadir, transcriptomefile)
sequencefile = os.path.join(datadir, sequencefile)

with open(transcriptomefilepath, mode='r') as infile:
    reader = csv.reader(infile)
    reader.next() # skip header line
    transcriptome = {rows[1]:int(float(rows[9].replace(',','.'))) for rows in reader} # (percentage of the) transcriptome
    #transcriptome = {rows[1]:10 for rows in reader} # constant transcriptome
print sum(transcriptome.values())

handle = open(sequencefile, "rU")
record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
handle.close()
orf_genomic_dict = {key:str(record_dict[key].seq).lower().replace('t', 'u') for key in record_dict}
#print orf_genomic_dict["YAL008W"] 

#pkl.dump(transcriptome, open("transcriptome.p", "wb")) # we now get the transcriptome from Premal Shah!
pkl.dump(orf_genomic_dict, open("orf_coding.p", "wb"))