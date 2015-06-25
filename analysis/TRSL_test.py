'''
Created on 23.04.2014

Runs the (non-specific) TRSL module standalone and writes the resulting timecourses to pickle files 

@author: MJS
'''
import os.path
import sys
import time
import pickle
import logging as log
import random as ran
from time import gmtime, strftime

sys.path.insert(0, './python_models/')
log.basicConfig(level=log.DEBUG, format='%(message)s', stream=sys.stdout)

import TRSL

def generate_timecourses_pickles(ribolist, n_mrnalist):
    for ribo in ribolist:
        for n_mrna in n_mrnalist:
            trsl = TRSL.TRSL()
            trsl.ribo_free = ribo
            trsl.n_mRNA = n_mrna
            trsl.gene_expressions = [ran.randint(1, 3795) for k in range(trsl.n_mRNA)]
            trsl.mRNAs = [TRSL.MRNA(index=gene) for gene in trsl.gene_expressions]
            trsl.solve_internal(0.0, 150.0, deltat=1.0)
            '''
            # Profiling:
            cProfile.run('trsl.solve(0.0, 60.0, deltat=1.0)', 'trsl_profile')
            p=pstats.Stats('trsl_profile')
            p.strip_dirs().sort_stats('cumulative').print_stats('TRSL')
            '''
            now = strftime("%Y%m%d_%H%M%S", gmtime())
            pickle_filename = now + "_" + str(ribo) + "_ribos_" + str(n_mrna) + "_mRNAs_timecourses_TRSL.pkl"
            pickle.dump({'trange':trsl.timerange, "timecourses":trsl.timecourses}, open(os.path.join(pickle_filename), "wb"))
            time.sleep(1.0)
            print

if __name__=="__main__":
    ribolist = [2, 20, 200, 2000]
    n_mrnalist = [6, 60, 600, 6000]
    generate_timecourses_pickles(ribolist, n_mrnalist)

    '''
    # test the ribosome ramp
    log.basicConfig(level=log.WARNING, format='%(message)s', stream=sys.stdout)
    nribo = 200000
    nmrna = 60
    trsl = TRSL.TRSL()
    trsl.ribo_free = nribo
    trsl.n_mRNA = nmrna
    trsl.gene_expressions = [ran.randint(1, 3795) for k in range(trsl.n_mRNA)]
    trsl.mRNAs = [TRSL.MRNA(index=gene) for gene in trsl.gene_expressions]
    trsl.solve_internal(0.0, 300.0, deltat=1.0)
    #for mRNA in trsl.mRNAs:
    #    print mRNA.ribosomes
    #print sum([len(mRNA.ribosomes.keys()) for mRNA in trsl.mRNAs])
    print
    now = strftime("%Y%m%d_%H%M%S", gmtime())
    pickle_filename = now + "_ribosomes_TRSL.pkl"
    pickle.dump([mRNA.ribosomes for mRNA in trsl.mRNAs], open(os.path.join(pickle_filename), "wb"))
    '''