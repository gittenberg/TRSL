'''
Script to simulate a grid of the following conditions:

- transcriptomes estimated at early G1 and S phases
- 100000, 200000 ad 30000 ribosomes

i.e. in total a grid of 2 x 3 = 6 conditions.

Created on 22.04.2016

@author: MJS
'''
print "TRSL_simulate_G1_S_transcriptome: starting"
import logging as log
import sys
import cPickle as pkl
import collections as col

from translation import MRNA_specific, TRSL_specific

log.basicConfig(level=log.DEBUG, format='%(message)s', stream=sys.stdout)

# load the two transcriptomes
phases = {0: 'early G1', 25: 'S'}
transcriptomes_dict = pkl.load((open('../parameters/transcriptome_time_dependent.p')))
transcriptomes_dict = {phases[key]: transcriptomes_dict[key] for key in transcriptomes_dict if key==0 or key==25}  # 0: early G1, 25: S phase

# load other data
exome = pkl.load(open("../parameters/orf_coding.p", "rb"))
init_rates = pkl.load(open("../parameters/init_rates_enhanced_median.p", "rb"))  # missing replaced by median

duration = 1200.0  # should be sufficent to saturate
ribonumbers = [100000, 200000, 300000]

for ribonumber in ribonumbers:
    for phase_ID in transcriptomes_dict:
        mRNAs = []
        description = '{} ribosomes, {} phase transcriptome, full exome, no decay, median-enhanced initiation rates according to Shah'.format(ribonumber, phase_ID)

        counter = 0
        genes = list(set(exome) & set(transcriptomes_dict[phase_ID]) & set(init_rates))
        print "found %s genes in common." % len(genes)

        for gene in genes:
            for transcript in range(transcriptomes_dict[phase_ID][gene]):
                mRNAs.append(MRNA_specific.mRNA_spec(index=counter, sequence=exome[gene], geneID=gene, ribosomes={}, init_rate=init_rates[gene]))  # do not just multiply the list
                counter += 1

        print "created transcriptome: {}.".format(description)

        tr = TRSL_specific.TRSL_spec(mRNAs, exome, decay_constants=None, nribo=ribonumber, proteome=col.Counter({}), detail=False)

        # tr._tRNA = col.Counter({i: TRSL_specific.tRNA_types[i]['abundancy'] for i in TRSL_specific.tRNA_types})
        tr._tRNA = col.Counter({i: int(TRSL_specific.tRNA_types[i]['abundancy'] * len(genes) / len(exome)) for i in TRSL_specific.tRNA_types})  # we do not let tRNA vary like the ribosomes
        tr._tRNA_free = col.Counter({i: int(tr._tRNA[i]) for i in TRSL_specific.tRNA_types})  # tRNA not bound to ribosomes
        tr._tRNA_bound = tr._tRNA - tr._tRNA_free  # tRNA bound to ribosomes

        # Run without profiling:
        # tr.solve_internal(0, duration, deltat=0.2)

        # Run with profiling:
        import cProfile
        cProfile.run('tr.solve_internal(0, '+str(duration)+', deltat=0.2)', 'trsl_profile')
        import pstats
        p = pstats.Stats('trsl_profile')
        p.strip_dirs().sort_stats('cumulative').print_stats()

        tr.dump_results(description)

