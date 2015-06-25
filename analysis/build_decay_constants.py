'''
Created on 28.01.2015

@author: martin
'''
import math
import os.path
import csv
import cPickle as pkl

datadir = "../data"
# http://www.pnas.org/content/103/35/13004.full
# half_lifes in minutes
half_lifes_file = r"Protein_half_lifes_SuppDataSet.csv"

half_lifes_filepath = os.path.join(datadir, half_lifes_file)
decay_constants = {}

with open(half_lifes_filepath, mode='r') as infile:
    reader = csv.reader(infile, delimiter=';')
    reader.next() # skip header line
    for rows in reader:
        try:
            decay_constants[rows[0]] = math.log(2.)/(float(rows[3])*60.) # in s^-1
        except:
            pass
            #print "found missing value, proceeding..."

#print decay_constants
print len(decay_constants)
pkl.dump(decay_constants, open("decay_constants.p", "wb"))
