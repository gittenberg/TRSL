"""
Script to simulate a cell cycle with varying transcriptomes

The cell cycle genes Sic1, Cln2 and Clb5 are turned on and off in different phases of the cell cycle

# This file is obsolete but might be reused for running another "dynamical transcriptome"
# TODO: Einschwingvorgang

__author__ = 'martin'
"""

import cPickle as pkl
import collections as col
from translation import MRNA_specific, TRSL_specific

transcriptomes_dict = pkl.load((open('../parameters/transcriptomes_cell_cycle.p')))

switch_times = sorted(transcriptomes_dict.keys())
switch_times.append(129 * 60)  # cell cycle takes 129 minutes

exome = pkl.load(open("../parameters/orf_coding.p", "rb"))
init_rates = pkl.load(open("../parameters/init_rates_plotkin.p", "rb"))
decay_constants = pkl.load(open("../parameters/decay_constants.p", "rb"))
transcriptome = transcriptomes_dict[0]
nribo = 200000

genes = list(set(exome) & set(transcriptome) & set(init_rates) & set(decay_constants))
print "{} genes found.".format(len(genes))

for start, stop in zip(switch_times[:-1], switch_times[1:]):
    print "simulating from {} to {}...".format(start, stop)

    transcriptome = transcriptomes_dict[start]
    mRNAs = []
    counter = 0
    description = 'polyphasic cell cycle from {} to {}, updated Shah transcriptome, full exome, no decay, updated initiation rates according to Shah'.format(start, stop)

    for gene in genes:
        for transcript in range(transcriptome[gene]):
            mRNAs.append(MRNA_specific.mRNA_spec(index=counter, sequence=exome[gene], geneID=gene, ribosomes={}, init_rate=init_rates[gene]))  # do not just multiply the list
            counter += 1

    print "created transcriptome at time {}.".format(start)

    tr = TRSL_specific.TRSL_spec(mRNAs, exome, decay_constants, nribo, detail=True)

    tr._tRNA = col.Counter({i: TRSL_specific.tRNA_types[i]['abundancy'] for i in TRSL_specific.tRNA_types})
    tr._tRNA_free = col.Counter({i: int(tr._tRNA[i]) for i in TRSL_specific.tRNA_types})  # tRNA not bound to ribosomes
    tr._tRNA_bound = tr._tRNA - tr._tRNA_free  # tRNA bound to ribosomes

    # TODO: carry over correct number of ribosomes

    # Profiling:
    import cProfile
    cProfile.run('tr.solve_internal('+str(start)+', '+str(stop)+', deltat=1.0)', 'trsl_profile')
    import pstats
    p=pstats.Stats('trsl_profile')
    p.strip_dirs().sort_stats('cumulative').print_stats()

    tr.dump_results(description)
