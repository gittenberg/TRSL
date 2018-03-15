'''
Script to simulate a grid of the following conditions:

- transcriptomes estimated at S phase and scaled down by a factor of 2
- different numbers of ribosomes
- this version with a constant transcriptome as per Ana's suggestion

to test the hypothesis of an optimal translational efficiency

Created on 15.03.2018

@author: MJS
'''
import logging as log
import sys
import cPickle as pkl
import collections as col
import numpy as np

from translation import MRNA_specific, TRSL_specific

print "TRSL_vary_ribosomes_tRNA_uniform_transcriptome: starting"
log.basicConfig(level=log.DEBUG, format='%(message)s', stream=sys.stdout)

# load the two transcriptomes
#######################################################################################################################
# setup 2: average transcriptome (60000 transcripts) and 50% of it
transcriptome_plotkin = pkl.load((open('../parameters/transcriptome_plotkin.p')))

# we now just modify the plotkin transcriptome to make it constant
ngenes = len(transcriptome_plotkin)  # 4839
ntranscripts = sum(transcriptome_plotkin.values())  # 60000

transcriptome_plotkin = {key: ntranscripts/ngenes for key in transcriptome_plotkin}  # yes integer division, always 12

transcriptome_scaled = {key: int(transcriptome_plotkin[key]*0.534) for key in transcriptome_plotkin}  # the 0.534 to compensate for the rounding down
transcriptomes_dict = {60000: transcriptome_plotkin, 30061: transcriptome_scaled}

for transcriptome in transcriptomes_dict:
    print "found {} transcripts in transcriptome {}...".format(sum(transcriptomes_dict[transcriptome].values()), transcriptome)

# load other data
# we now just modify the exome and give it a standard transcript
exome = pkl.load(open("../parameters/orf_coding.p", "rb"))
standard_transcript = 'aug' + (1251 - 6) / 3 * 'ggu' + 'uaa'  # ggu because it has a lot of tRNAs
exome = {key: standard_transcript for key in exome}

# we now just modify the plotkin init rates and replace them by their median
init_rates = pkl.load(open("../parameters/init_rates_enhanced_median.p", "rb"))  # missing replaced by median
median_init_rate = np.median(init_rates.values())  # 1.29e-06
init_rates = {key: median_init_rate for key in init_rates}

duration = 1200.0  # should be sufficent to saturate
ribonumbers = range(50000, 550000, 50000)

for ribonumber in ribonumbers:
    scaling_factor = ribonumber * 1.0 / 200000
    for transcriptome_ID in transcriptomes_dict:
        mRNAs = []
        description = '{} ribosomes, {} phase transcriptome, varying tRNAs, constant exome, no decay, median constant initiation rates, deltat=0.2s'.format(ribonumber, transcriptome_ID)

        counter = 0
        genes = list(set(exome) & set(transcriptomes_dict[transcriptome_ID]) & set(init_rates))
        print "found %s genes in common." % len(genes)

        for gene in genes:
            for transcript in range(transcriptomes_dict[transcriptome_ID][gene]):
                mRNAs.append(MRNA_specific.mRNA_spec(index=counter, sequence=exome[gene], geneID=gene, ribosomes={}, init_rate=init_rates[gene]))  # do not just multiply the list
                counter += 1

        print "created transcriptome: {}.".format(description)

        tr = TRSL_specific.TRSL_spec(mRNAs, exome, decay_constants=None, nribo=ribonumber, proteome=col.Counter({}), detail=False)

        # tr._tRNA = col.Counter({i: TRSL_specific.tRNA_types[i]['abundancy'] for i in TRSL_specific.tRNA_types})
        tr._tRNA = col.Counter({i: int(TRSL_specific.tRNA_types[i]['abundancy'] * scaling_factor * len(genes) / len(exome)) for i in TRSL_specific.tRNA_types})  # this time we do let tRNA vary like the ribosomes
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
