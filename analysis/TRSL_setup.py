'''
Created on 29.01.2015

@author: MJS
'''
print "TRSL_setup: starting"
import logging as log
import sys
import cPickle as pkl
import collections as col

from translation import MRNA_specific, TRSL_specific

log.basicConfig(level=log.DEBUG, format='%(message)s', stream=sys.stdout)

examplesequence_1 = "ggg uuu uca uca uuu gag gac gau gua aaa ggg uaa".replace(' ', '')
examplesequence_2 = "aug aaa cug ccc gag ggg uuu uca uca uuu gag gac aaa cug ccc gag ggg uuu uca uca uuu gag gac aaa cug ccc gag ggg uuu uca uca uaa".replace(' ', '')

# configuration dictionary
conf = {1: {
    'exome': {1: examplesequence_1, 2: examplesequence_2},
    'transcriptome': {1: 2, 2: 1},
    'init_rates': {1: 0.1, 2: 0.1},
    'description': 'test configuration with 2 genes in 3 transcripts'
}, 2: {
    'exome': pkl.load(open("../parameters/orf_coding.p", "rb")),
    'transcriptome': pkl.load(open("../parameters/transcriptome_plotkin_20000.p", "rb")),
    'init_rates': pkl.load(open("../parameters/init_rates_plotkin.p", "rb")),
    'description': '20000 transcriptome, full exome, no decay, Plotkin initiation probabilities'
}, 3: {
    'exome': pkl.load(open("../parameters/orf_coding.p", "rb")),
    'transcriptome': pkl.load(open("../parameters/transcriptome_plotkin.p", "rb")),
    'init_rates': {gene: 8.2e-07 for gene in pkl.load(open("../parameters/init_rates_plotkin.p", "rb"))},
    # average initiation rate
    'description': 'full transcriptome and exome, no decay, constant initiation rates'
}, 4: {
    'exome': pkl.load(open("../parameters/orf_coding.p", "rb")),
    'transcriptome': pkl.load(open("../parameters/transcriptome_plotkin.p", "rb")),
    'init_rates': pkl.load(open("../parameters/init_rates_plotkin.p", "rb")),
    'description': 'full transcriptome and exome, no decay, specific best estimate initiation rates according to Plotkin'
}, 5: {
    'exome': pkl.load(open("../parameters/orf_coding.p", "rb")),
    'transcriptome': pkl.load(open("../parameters/transcriptome_plotkin.p", "rb")),
    'init_rates': pkl.load(open("../parameters/init_rates_plotkin.p", "rb")),
    'decay_constants': pkl.load(open("../parameters/decay_constants.p", "rb")),
    'description': 'full transcriptome and exome, with decay, specific best estimate initiation rates according to Plotkin'
}, 6: {
    'exome': pkl.load(open("../parameters/orf_coding.p", "rb")),
    'transcriptome': pkl.load(open("../parameters/transcriptome.p", "rb")),
    'init_rates': pkl.load(open("../parameters/init_rates_stansfield.p", "rb")),
    'decay_constants': pkl.load(open("../parameters/decay_constants.p", "rb")),
    'description': 'full transcriptome and exome, specific best estimate initiation rates according to Stansfield'
}, 7: {
    'exome': pkl.load(open("../parameters/orf_coding.p", "rb")),
    'transcriptome': pkl.load(open("../parameters/transcriptome_plotkin.p", "rb")),
    'init_rates': pkl.load(open("../parameters/init_rates_plotkin_old_1.p", "rb")),
    'description': 'full transcriptome and exome, no decay, old (buggy) initiation rates according to Plotkin'
}, 8: {
    # base case run with reduced number of genes (limited by availability of initiation rates)
    'exome': pkl.load(open("../parameters/orf_coding.p", "rb")),
    'transcriptome': pkl.load(open("../parameters/transcriptome_shah.p", "rb")),
    'init_rates': pkl.load(open("../parameters/init_rates_plotkin.p", "rb")),
    'description': 'updated Shah transcriptome, full exome, no decay, updated initiation rates according to Shah'
}, 9: {
    'exome': pkl.load(open("../parameters/orf_coding.p", "rb")),
    'transcriptome': pkl.load(open("../parameters/transcriptome_shah.p", "rb")),
    'init_rates': pkl.load(open("../parameters/init_rates_plotkin.p", "rb")),
    'decay_constants': pkl.load(open("../parameters/decay_constants.p", "rb")),
    'description': 'updated Shah transcriptome, full exome, with decay, updated initiation rates according to Shah'
}, 10: {
    'exome': pkl.load(open("../parameters/orf_coding.p", "rb")),
    'transcriptome': pkl.load(open("../parameters/transcriptome_shah.p", "rb")),
    'init_rates': pkl.load(open("../parameters/init_rates_plotkin.p", "rb")),
    'decay_constants': pkl.load(open("../parameters/decay_constants.p", "rb")),
    'description': 'updated Shah transcriptome, full exome, with decay, rare tRNAs reduced 90%, updated initiation rates according to Shah',
    'tRNA': col.Counter(
            {i: TRSL_specific.tRNA_types[i]['abundancy'] if TRSL_specific.tRNA_types[i]['abundancy'] > 11070 else 1107
             for i in TRSL_specific.tRNA_types})
}, 11: {
    'exome': pkl.load(open("../parameters/orf_coding.p", "rb")),
    'transcriptome': pkl.load(open("../parameters/transcriptome_shah.p", "rb")),
    'init_rates': pkl.load(open("../parameters/init_rates_plotkin_old_1.p", "rb")),
    #'decay_constants': pkl.load(open("../parameters/decay_constants.p", "rb")),
    'description': 'updated Shah transcriptome, full exome, without decay, OLD initiation rates according to Shah',
}}

if __name__ == "__main__":
    log.basicConfig(level=log.DEBUG, format='%(message)s', stream=sys.stdout)

    for i in [8]:  # set configuration_id
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
                    # print "abundancies and initiation rates available for gene:", gene
                    for instance in range(conf[i]['transcriptome'][gene]):
                        mRNAs.append(MRNA_specific.mRNA_spec(index=counter, sequence=conf[i]['exome'][gene], geneID=gene, ribosomes={}, init_rate=conf[i]['init_rates'][gene]))  # do not just multiply the list
                        counter += 1
            else:
                if gene in conf[i]['transcriptome']:
                    print "abundancies but no initiation rate available for gene:", gene
                    for instance in range(conf[i]['transcriptome'][gene]):
                        mRNAs.append(MRNA_specific.mRNA_spec(index=counter, sequence=conf[i]['exome'][gene], geneID=gene, ribosomes={}))  # do not just multiply the list
                        counter += 1
        print "built gene library, next: run TRSL_spec."

        description = conf[i]['description']
        print description

        duration = 1800.0

        tr = TRSL_specific.TRSL_spec(mRNAs, conf[i]['exome'], conf[i]['decay_constants'], nribo=200000, detail=True)

        # overwrite tRNA:
        if 'tRNA' not in conf[i]:
            tr._tRNA = col.Counter({i: TRSL_specific.tRNA_types[i]['abundancy'] for i in TRSL_specific.tRNA_types}) # full tRNA
        else:
            tr._tRNA = conf[i]['tRNA']
        # print tr._tRNA
        tr._tRNA_free = col.Counter({i: int(tr._tRNA[i]) for i in TRSL_specific.tRNA_types})  # tRNA not bound to ribosomes
        tr._tRNA_bound = tr._tRNA - tr._tRNA_free  # tRNA bound to ribosomes

        print "found {} mRNAs".format(len(mRNAs))
        print "found {} genes".format(len(conf[i]['exome']))
        print "found {} tRNA molecules".format(sum(tr._tRNA.values()))
        print "solving..."

        tr.solve_internal(0.0, duration, deltat=0.2)

        '''
        # Profiling:
        import cProfile
        cProfile.run('tr.solve_internal(0.0, '+str(duration)+', deltat=0.1)', 'trsl_profile')
        import pstats
        p = pstats.Stats('trsl_profile')
        p.strip_dirs().sort_stats('cumulative').print_stats()
        #tr.inspect()

        tr.dump_results(description)

        # write last polysomes to shelve database
        tr.write_last_polysome(description)
        '''
    print "done."
