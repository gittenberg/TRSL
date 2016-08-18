'''
Created on 22.02.2015

@author: MJS
'''
import os.path
import csv
import cPickle as pkl

datadir = "../data"
# http://www.proteomaps.net/data_sets/sce_Nagaraj/sce_Nagaraj.csv
annotations_file = r"../data/proteomaps_rich_medium_1_sce_Nagaraj.sbt"
annotations_filepath = os.path.join(datadir, annotations_file)

annotations = {}
proteome = {}
with open(annotations_filepath, mode='r') as infile:
    reader = csv.reader(infile)
    reader.next() # skip header line
    reader.next() # skip second line
    for rows in reader:
        row = rows[0].split('\t')
        abundance = float(row[1].replace(',','.'))
        if abundance>1e-20:
            proteome[row[0]] = abundance
        annotations[row[0]] = row[8]

print annotations
print len(annotations)
print proteome
print len(proteome)
pkl.dump(annotations, open("../parameters/annotations_liebermeister.p", "wb"))
pkl.dump(proteome, open("../parameters/prot_nagaraj.p", "wb"))
