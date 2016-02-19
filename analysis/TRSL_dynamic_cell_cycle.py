"""
Script to simulate a cell cycle with varying transcriptomes

__author__ = 'martin'
"""

import cPickle as pkl
import collections as col
import numpy as np
from translation import MRNA_specific, TRSL_specific

transcriptomes_dict = pkl.load((open('../parameters/transcriptome_time_dependent.p')))

# when are new transcriptomes loaded (conversion from minutes to seconds):
switch_times = [key * 60 for key in sorted(transcriptomes_dict.keys())]

# load other data
exome = pkl.load(open("../parameters/orf_coding.p", "rb"))
# init_rates = pkl.load(open("../parameters/init_rates_plotkin.p", "rb"))
init_rates = pkl.load(open("../parameters/init_rates_enhanced_median.p", "rb"))  # missing replaced by median
decay_constants = pkl.load(open("../parameters/decay_constants.p", "rb"))
# load initial transcriptome
transcriptome = transcriptomes_dict[0]

# find common data set
# genes = list(set(exome) & set(transcriptome) & set(init_rates) & set(decay_constants)) # with decay
genes = list(set(exome) & set(transcriptome) & set(init_rates)) # without decay
print "{} genes found.".format(len(genes))

# to create a growing number of ribosomes and tRNAs
nribo_start = 200000 * len(genes) / len(exome)  # scaled to make ribosome count more realistic
growth_factor_range = np.arange(1., 1.5, (1.5-1.0)/len(switch_times))  # by how much they grow in each interval

# run simulation
# Einschwingvorgang: reichen 900 s?
burnin = 900
for start, stop, growth_factor in zip(switch_times[:-1], switch_times[1:], growth_factor_range):
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

    tr = TRSL_specific.TRSL_spec(mRNAs, exome, decay_constants, int(nribo_start*growth_factor), proteome=col.Counter({}), detail=True)

    # tr._tRNA = col.Counter({i: TRSL_specific.tRNA_types[i]['abundancy'] for i in TRSL_specific.tRNA_types})
    tr._tRNA = col.Counter({i: int(TRSL_specific.tRNA_types[i]['abundancy'] * len(genes) / len(exome) * growth_factor) for i in TRSL_specific.tRNA_types})
    tr._tRNA_free = col.Counter({i: int(tr._tRNA[i]) for i in TRSL_specific.tRNA_types})  # tRNA not bound to ribosomes
    tr._tRNA_bound = tr._tRNA - tr._tRNA_free  # tRNA bound to ribosomes

    # Run without profiling:
    tr.solve_internal(start, stop+burnin, deltat=1.0)
    '''
    # Run with profiling:
    import cProfile
    cProfile.run('tr.solve_internal('+str(start)+', '+str(stop+900)+', deltat=1.0)', 'trsl_profile') # 900 s burn-in
    import pstats
    p=pstats.Stats('trsl_profile')
    p.strip_dirs().sort_stats('cumulative').print_stats()
    '''
    tr.dump_results(description)
