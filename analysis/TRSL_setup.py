'''
Created on 29.01.2015

@author: MJS
'''
print "TRSL_setup: starting"
import logging as log
import sys
import cPickle as pkl
import collections as col
import cProfile
import pstats

#import MRNA_specific  # @UnresolvedImport
#import TRSL_specific

from translation import MRNA_specific, TRSL_specific

log.basicConfig(level=log.DEBUG, format='%(message)s', stream=sys.stdout)

examplesequence_1 = "ggg uuu uca uca uuu gag gac gau gua aaa ggg uaa".replace(' ', '')
examplesequence_2 = "aug aaa cug ccc gag ggg uuu uca uca uuu gag gac aaa cug ccc gag ggg uuu uca uca uuu gag gac aaa cug ccc gag ggg uuu uca uca uaa".replace(' ', '')

# configuration dictionary
conf = {}
conf[1] = {
           'exome': {1: examplesequence_1, 2: examplesequence_2},
           'transcriptome': {1: 1, 2: 2},
           'init_rates': {1: 0.1, 2: 0.1},
           'description': 'test configuration with 2 genes in 3 transcripts'
           }
conf[2] = {
           'exome': pkl.load(open("../parameters/orf_coding.p", "rb" )),
           'transcriptome': pkl.load(open("../parameters/transcriptome.p", "rb" )),
           'init_rates': pkl.load(open("../parameters/init_rates_stansfield.p", "rb" )),
           'description': 'full transcriptome and exome, specific best estimate initiation rates according to Stansfield'
           }
conf[3] = {
           'exome': pkl.load(open("../parameters/orf_coding.p", "rb" )),
           'transcriptome': pkl.load(open("../parameters/transcriptome.p", "rb" )),
           'init_rates': {gene: 8.2e-07 for gene in pkl.load(open("../parameters/init_rates_plotkin.p", "rb" ))}, # average initiation rate
           'description': 'full transcriptome and exome, no decay, constant initiation rates'
           }
conf[4] = {
           'exome': pkl.load(open("../parameters/orf_coding.p", "rb" )),
           'transcriptome': pkl.load(open("../parameters/transcriptome.p", "rb" )),
           'init_rates': pkl.load(open("../parameters/init_rates_plotkin.p", "rb" )),
           'description': 'full transcriptome and exome, no decay, specific best estimate initiation rates according to Plotkin'
           }
conf[5] = {
           'exome': pkl.load(open("../parameters/orf_coding.p", "rb" )),
           'transcriptome': pkl.load(open("../parameters/transcriptome.p", "rb" )),
           'init_rates': pkl.load(open("../parameters/init_rates_plotkin.p", "rb" )),
           'decay_constants': pkl.load(open("../parameters/decay_constants.p", "rb" )),
           'description': 'full transcriptome and exome, with decay, specific best estimate initiation rates according to Plotkin'
           }
conf[6] = {
           'exome': pkl.load(open("../parameters/orf_coding.p", "rb" )),
           'transcriptome': pkl.load(open("../parameters/transcriptome.p", "rb" )),
           'init_rates': pkl.load(open("../parameters/init_rates_stansfield.p", "rb" )),
           'decay_constants': pkl.load(open("../parameters/decay_constants.p", "rb" )),
           'description': 'full transcriptome and exome, specific best estimate initiation rates according to Stansfield'
           }

for i in [5]: # set configuration_id
    if 'decay_constants' in conf[i]:
        genes = list(set(conf[i]['exome']) & set(conf[i]['transcriptome']) & set(conf[i]['init_rates']) & set(conf[i]['decay_constants']))
    else:
        genes = list(set(conf[i]['exome']) & set(conf[i]['transcriptome']) & set(conf[i]['init_rates']))
        conf[i]['decay_constants'] = None
    print "found %s genes in common." % len(genes)
    mRNAs = []
    counter = 0
    for gene in genes:
        if conf[i]['init_rates']:
            if gene in conf[i]['transcriptome'] and gene in conf[i]['init_rates']:
                print "abundancies and initiation rates available for gene:", gene
                for instance in range(conf[i]['transcriptome'][gene]):
                    mRNAs.append(MRNA_specific.mRNA_spec(index=counter, sequence=conf[i]['exome'][gene], geneID=gene, ribosomes={3: None, 9: None}, init_rate=conf[i]['init_rates'][gene])) # do not just multiply the list
                    counter += 1
        else:
            if gene in conf[i]['transcriptome']:
                print "abundancies but no initiation rate available for gene:", gene
                for instance in range(conf[i]['transcriptome'][gene]):
                    mRNAs.append(MRNA_specific.mRNA_spec(index=counter, sequence=conf[i]['exome'][gene], geneID=gene, ribosomes={3:None, 9:None})) # do not just multiply the list
                    counter += 1
    print "built gene library, next: run TRSL_spec."
    
    print conf[i]['description']
    
    duration = 49.0
    
    tr = TRSL_specific.TRSL_spec(mRNAs, conf[i]['exome'], conf[i]['decay_constants'], proteome={"YKR017C":0})
    #tr.tRNA = col.Counter({i:TRSL_specific.tRNA_types[i]['abundancy']*2 for i in TRSL_specific.tRNA_types}) # double tRNA inventory to prevent stalling
    tr.tRNA = col.Counter({i:TRSL_specific.tRNA_types[i]['abundancy'] for i in TRSL_specific.tRNA_types}) # do not double tRNA inventory to induce stalling
    tr.tRNA_free = col.Counter({i:int(tr.tRNA[i]) for i in TRSL_specific.tRNA_types})     # tRNA not bound to ribosomes
    tr.tRNA_bound = tr.tRNA - tr.tRNA_free                                # tRNA bound to ribosomes
    tr.solve_internal(0.0, duration, deltat=1.0)
    
    '''
    # Profiling:
    cProfile.run('tr.solve_internal(0.0, '+str(duration)+', deltat=1.0)', 'trsl_profile')
    p=pstats.Stats('trsl_profile')
    p.strip_dirs().sort_stats('cumulative').print_stats()
    #tr.inspect()
    '''

    results = {}
    results['proteome'] = tr.proteins
    results['timerange'] = tr.timerange
    results['timecourses'] = tr.timecourses
    results["description"] = conf[i]['description']
    import time; now = time.strftime("%Y%m%d_%H%M", time.gmtime())
    results["time_stamp"] = now
    results["n_ribosomes"] = tr.ribo_bound + tr.ribo_free
    results["n_tRNA"] = sum(tr.tRNA.values())
    results["duration"] = duration
    pkl.dump(results, open("../results/results_"+now+"_"+str(int(duration)).zfill(4)+"s.p", "wb"))
    print conf[i]['description']

print "done."
