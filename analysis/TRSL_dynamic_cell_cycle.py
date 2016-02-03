"""
Script to simulate a cell cycle with varying transcriptomes

# This file is obsolete but might be reused for running another "dynamical transcriptome"
# TODO: Einschwingvorgang

__author__ = 'martin'
"""

import cPickle as pkl
import collections as col
from translation import MRNA_specific, TRSL_specific

transcriptomes_dict = pkl.load((open('../parameters/transcriptome_time_dependent.p')))

# when are new transcriptomes loaded:
switch_times = [key * 60 for key in sorted(transcriptomes_dict.keys())]

# load other data and initial transcriptome
exome = pkl.load(open("../parameters/orf_coding.p", "rb"))
init_rates = pkl.load(open("../parameters/init_rates_plotkin.p", "rb"))
decay_constants = pkl.load(open("../parameters/decay_constants.p", "rb"))
transcriptome = transcriptomes_dict[0]

# find common data set
#genes = list(set(exome) & set(transcriptome) & set(init_rates) & set(decay_constants)) # with decay # TODO: only 3000-ish genes: adjust
genes = list(set(exome) & set(transcriptome) & set(init_rates)) # without decay # TODO: only 4682-ish genes: adjust
print "{} genes found.".format(len(genes))

# create a growing number of ribosomes
nribo_start = 200000 * len(genes) / len(exome) # scaled to make ribosome count more realistic
nsribo = range(nribo_start, int(1.5*nribo_start), int((1.5-1.0)*nribo_start/len(switch_times)))

# run simulation
# Einschwingvorgang: reichen 900 s?
burnin = 900
for start, stop, nribo in zip(switch_times[:-1], switch_times[1:], nsribo):
    print "simulating from {} to {}...".format(start, stop)

    transcriptome = transcriptomes_dict[start/60]
    mRNAs = []
    counter = 0
    description = 'polyphasic cell cycle from {} to {}, updated Shah transcriptome, full exome, no decay, updated initiation rates according to Shah'.format(start, stop)

    for gene in genes:
        for transcript in range(transcriptome[gene]):
            mRNAs.append(MRNA_specific.mRNA_spec(index=counter, sequence=exome[gene], geneID=gene, ribosomes={}, init_rate=init_rates[gene]))  # do not just multiply the list
            counter += 1

    print "created transcriptome at time {}.".format(start)

    tr = TRSL_specific.TRSL_spec(mRNAs, exome, decay_constants, nribo, detail=True)

    #tr._tRNA = col.Counter({i: TRSL_specific.tRNA_types[i]['abundancy'] for i in TRSL_specific.tRNA_types})
    tr._tRNA = col.Counter({i: TRSL_specific.tRNA_types[i]['abundancy'] * len(genes) / len(exome) for i in TRSL_specific.tRNA_types})
    tr._tRNA_free = col.Counter({i: int(tr._tRNA[i]) for i in TRSL_specific.tRNA_types})  # tRNA not bound to ribosomes
    tr._tRNA_bound = tr._tRNA - tr._tRNA_free  # tRNA bound to ribosomes

    # Run without profiling:
    tr.solve_internal(start, stop+burnin, deltat=1.0)
    '''
    # Run with profiling:
    import cProfile
    cProfile.run('tr.solve_internal('+str(start)+', '+str(stop+300)+', deltat=1.0)', 'trsl_profile') # 300 s burn-in # TODO: needs to be removed after simulation
    import pstats
    p=pstats.Stats('trsl_profile')
    p.strip_dirs().sort_stats('cumulative').print_stats()
    '''
    tr.dump_results(description)
    tr.timecourses = {} # not sure if works or necessary
