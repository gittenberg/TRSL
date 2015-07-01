__author__ = 'Max Floettmann'

import logging as log
import sys
import cPickle as pkl
import collections as col

from translation import MRNA_specific, TRSL_specific

log.basicConfig(level=log.DEBUG, format='%(message)s', stream=sys.stdout)

examplesequence_1 = "ggg uuu uca uca uuu gag gac gau gua aaa ggg uaa".replace(' ', '')
examplesequence_2 = "aug aaa cug ccc gag ggg uuu uca uca uuu gag gac aaa cug ccc gag ggg uuu uca uca uuu gag gac aaa cug ccc gag ggg uuu uca uca uaa".replace(' ', '')


conf = {}
conf[1] = dict(exome={1: examplesequence_1, 2: examplesequence_2},
               transcriptome={1: 1, 2: 2},
               init_rates={1: 0.1, 2: 0.1},
               description='test configuration with 2 genes in 3 transcripts')

conf[2] = dict(exome=pkl.load(open("../parameters/orf_coding.p", "rb")),
               transcriptome=pkl.load(open("../parameters/transcriptome.p", "rb")),
               init_rates={gene: 8.2e-07 for gene in pkl.load(open("../parameters/init_rates_plotkin.p", "rb"))},
               description='full transcriptome and exome, no decay, constant initiation rates')

for i in [2]:  # set configuration_id
    if 'decay_constants' in conf[i]:
        genes = list(set(conf[i]['exome']) & set(conf[i]['transcriptome']) & set(conf[i]['init_rates']) & set(
            conf[i]['decay_constants']))
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
                    mRNAs.append(MRNA_specific.mRNA_spec(index=counter,
                                                         sequence=conf[i]['exome'][gene],
                                                         geneID=gene,
                                                         ribosomes={},
                                                         init_rate=conf[i]['init_rates'][gene]))  # do not just multiply the list
                    counter += 1
        else:
            if gene in conf[i]['transcriptome']:
                print "abundancies but no initiation rate available for gene:", gene
                for instance in range(conf[i]['transcriptome'][gene]):
                    mRNAs.append(MRNA_specific.mRNA_spec(index=counter, sequence=conf[i]['exome'][gene], geneID=gene,
                                                         ribosomes={}))  # do not just multiply the list
                    counter += 1
    print "built gene library, next: run TRSL_spec."

    description = conf[i]['description']
    print description

    duration = 5.0

    tr = TRSL_specific.TRSL_spec(mRNAs, conf[i]['exome'], conf[i]['decay_constants'])

    proteins, mRNAs = tr.solve_internal(0.0, duration, deltat=1.0)

    for m in mRNAs:
        m.init_rate = 0
    tr.init_rate = 0

    tr.mRNAs = mRNAs
    proteins, mRNAs = tr.solve_internal(0.0, 5, deltat=1.0)

    '''
    # Profiling:
    cProfile.run('tr.solve_internal(0.0, '+str(duration)+', deltat=1.0)', 'trsl_profile')
    p=pstats.Stats('trsl_profile')
    p.strip_dirs().sort_stats('cumulative').print_stats()
    #tr.inspect()
    '''

    #tr.dump_results(description)

print "done."