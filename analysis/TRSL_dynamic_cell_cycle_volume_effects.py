"""
Script to simulate a cell cycle with varying transcriptomes
Includes effects of
- ribosome growth
- tRNA growth
- mRNA variation (time course)
- mRNA dilution (tRNA influence on elongation rate, ribosomes on initiation rate)

This script overrides TRSL_dynamic_cell_cycle.py

__author__ = 'martin'
"""

import cPickle as pkl
import numpy as np
import collections as col
from translation import MRNA_specific, TRSL_specific

transcriptomes_dict = pkl.load((open('../../parameters/transcriptome_time_dependent.p')))

# when are new transcriptomes loaded (conversion from minutes to seconds):
switch_times = [key * 60 for key in sorted(transcriptomes_dict.keys())]

# load exome
exome = pkl.load(open("../parameters/orf_coding.p", "rb"))
# load initiation rates
# init_rates = pkl.load(open("../parameters/init_rates_enhanced_median.p", "rb"))  # missing replaced by median
init_rates = pkl.load(open("../parameters/init_rates_plotkin.p", "rb"))
# load initial transcriptome
transcriptome = transcriptomes_dict[0]

total_growth = 1.7  # by which factor the cell grows during one cell cycle
# time-resolved growth factor at the beginnings of the intervals
growth_factor_range = np.arange(1.0, total_growth, (total_growth - 1.0) / len(switch_times))

# find common data set
genes = list(set(exome) & set(transcriptome) & set(init_rates))
print "TRSL_dynamic_cell_cycle_volume_effects: {} genes found.".format(len(genes))

# to create a growing number of ribosomes and tRNAs
nribo_start = 200000  # Needs no scaling as this is taken care of by transcriptome

# run simulation
# Einschwingvorgang: 1800 s
burnin = 1800
for start, stop, growth_factor in zip(switch_times[:-1], switch_times[1:], growth_factor_range):
    print "TRSL_dynamic_cell_cycle_volume_effects: simulating from {} to {}...".format(start, stop)

    transcriptome = transcriptomes_dict[start / 60]
    mRNAs = []
    counter = 0
    description = 'volume-adjusted polyphasic cell cycle from {} to {}, Teufel transcriptome, full exome, no decay, with ribo growth factor, updated initiation rates according to Shah'.format(start, stop)

    for gene in genes:
        for transcript in range(transcriptome[gene]):
            mRNAs.append(MRNA_specific.mRNA_spec(index=counter, sequence=exome[gene], geneID=gene, ribosomes={},
                                                 init_rate=init_rates[gene]))  # do not just multiply the list
            counter += 1

    print "TRSL_dynamic_cell_cycle_volume_effects: created transcriptome at time {}.".format(start)

    tr = TRSL_specific.TRSL_spec(mRNAs, exome, None, int(nribo_start * growth_factor), proteome=col.Counter({}), detail=True)

    tr._tRNA = col.Counter({i: int(TRSL_specific.tRNA_types[i]['abundancy'] * len(genes) / len(exome) * growth_factor) for i in TRSL_specific.tRNA_types})
    tr._tRNA_free = col.Counter({i: int(tr._tRNA[i]) for i in TRSL_specific.tRNA_types})  # tRNA not bound to ribosomes
    tr._tRNA_bound = tr._tRNA - tr._tRNA_free  # tRNA bound to ribosomes
    print "TRSL_dynamic_cell_cycle_volume_effects: corrected tRNA counts for growth effect"

    # reduce initiation rates of every mRNA by factor inversely proportional to volume
    for mRNA in tr.mRNAs:
        mRNA.init_rate /= growth_factor
    print "TRSL_dynamic_cell_cycle_volume_effects: corrected initation rates for growth effect"
    # reduce elongation rate by factor inversely proportional to volume
    tr.elong_rate /= growth_factor

    # Run without profiling:
    tr.solve_internal(start, stop + burnin, deltat=1.0)
    '''
    # Run with profiling:
    import cProfile
    cProfile.run('tr.solve_internal('+str(start)+', '+str(stop+900)+', deltat=1.0)', 'trsl_profile') # 900 s burn-in
    import pstats
    p=pstats.Stats('trsl_profile')
    p.strip_dirs().sort_stats('cumulative').print_stats()
    '''
    tr.dump_results(description)

print "done."
