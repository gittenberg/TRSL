'''
Created on 30.10.2014

@author: martin
'''
import TRSL_specific
import MRNA_specific
import logging as log
import sys
import cProfile
import pstats
import pickle
import random
import os.path

#log.basicConfig(level=log.DEBUG, format='%(message)s', stream=sys.stdout)
duration = 100.0

'''
examplesequence = "gug ucu ugc uau ucg agu gca ucc aga uca uac aag ugu aua cuu aug cgu gua gac cug cga cuc cua auu gau ggu guc \
                   guu gaa uuu ugg cgc auc cgg gcc ggg gga ggc gag acg uug uua ccg uuc acc aca aaa gcg aac ccu cag gcu agc cau \
                   aau agg acu cac caa cca ccc uaa".replace(' ', '')
gene_library = {1: examplesequence}
mRNAs = [TRSL_specific.mRNA_spec(index=0, sequence=examplesequence, geneID=1)]*60
'''

gene_library = pickle.load(open("gene_library.pkl", "r")) 

def build_transcriptome(transcriptome_size, gene_library):
    print "building transcriptome..."
    #print len(gene_library)
    transcribed_genes = [random.choice(gene_library.items()) for n in range(transcriptome_size)]
    transcribed_genes = [(name, sequence.lower().replace('t', 'u')) for name, sequence in transcribed_genes]
    mRNAs = [MRNA_specific.mRNA_spec(index=n, sequence=transcribed_genes[n][1], geneID=transcribed_genes[n][0]) for n in range(transcriptome_size)]
    pickle.dump(mRNAs, open(os.path.join("mRNA_"+str(transcriptome_size)+".pkl"), "wb"))
    return mRNAs

mRNAs = build_transcriptome(600, gene_library)
#mRNAs = pickle.load(open("mRNA_60000.pkl", "r")) 

tr = TRSL_specific.TRSL_spec(mRNAs, gene_library)
#tr.solve_internal(0.0, duration, deltat=1.0)

# Profiling:
cProfile.run('tr.solve_internal(0.0, duration, deltat=1.0)', 'trsl_profile')
p = pstats.Stats('trsl_profile')
p.strip_dirs().sort_stats('cumulative').print_stats(25)
