'''
Script to simulate a grid of the following conditions:

- transcriptomes estimated at S phase and scaled down by a factor of 2
- different numbers of ribosomes
- this version keeps tRNA constant

to test the hypothesis of an optimal translational efficiency

Created on 27.04.2016

Updated to cover wider range of ribosome numbers.

@author: MJS
'''
print "TRSL_vary_ribosomes: starting"
import logging as log
import sys
import cPickle as pkl
import collections as col

from translation import MRNA_specific, TRSL_specific

log.basicConfig(level=log.DEBUG, format='%(message)s', stream=sys.stdout)

# load the two transcriptomes
#######################################################################################################################
# setup 1: S phase transcriptome (ca. 72000 transcripts) and 50% of it
'''
transcriptomes_dict = pkl.load((open('../parameters/transcriptome_time_dependent.p')))
transcriptome_scaled = {key: int(round(transcriptomes_dict[25][key]*0.5, ndigits=0)) for key in transcriptomes_dict[25]}  # S phase is at 25 minutes
transcriptomes_dict = {'S': transcriptomes_dict[25], 'S_scaled': transcriptome_scaled}
'''
#######################################################################################################################
# setup 2: average transcriptome (60000 transcripts) and 50% of it
transcriptome_plotkin = pkl.load((open('../parameters/transcriptome_plotkin.p')))
transcriptome_scaled = {key: int(transcriptome_plotkin[key]*0.534) for key in transcriptome_plotkin}  # the 0.534 to compensate for the rounding down
transcriptomes_dict = {60000: transcriptome_plotkin, 30061: transcriptome_scaled}

for transcriptome in transcriptomes_dict:
    print "found {} transcripts in transcriptome {}...".format(sum(transcriptomes_dict[transcriptome].values()), transcriptome)

# load other data
exome = pkl.load(open("../parameters/orf_coding.p", "rb"))
init_rates = pkl.load(open("../parameters/init_rates_enhanced_median.p", "rb"))  # missing replaced by median
# init_rates = pkl.load(open("../parameters/init_rates_plotkin.p", "rb"))
#  TODO: check if this makes a difference; it should not because the number of common genes is the same both ways

duration = 1200.0  # 1200. should be sufficent to saturate
# ribonumbers = range(50000, 550000, 50000)
ribonumbers = [1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000]

for ribonumber in ribonumbers:

    for transcriptome_ID in transcriptomes_dict:
        mRNAs = []
        deltat = 0.05
        description = '{} ribosomes, {} phase transcriptome, constant tRNAs, full exome, median-enhanced initiation rates according to Shah, deltat={}s'.format(ribonumber, transcriptome_ID, deltat)

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
        tr._tRNA = col.Counter({i: int(TRSL_specific.tRNA_types[i]['abundancy'] * 1.0 * len(genes) / len(exome)) for i in TRSL_specific.tRNA_types})  # we do not let tRNA vary like the ribosomes
        tr._tRNA_free = col.Counter({i: int(tr._tRNA[i]) for i in TRSL_specific.tRNA_types})  # tRNA not bound to ribosomes
        tr._tRNA_bound = tr._tRNA - tr._tRNA_free  # tRNA bound to ribosomes

        # Run without profiling:
        # tr.solve_internal(0, duration, deltat=deltat)

        # Run with profiling:
        import cProfile
        filename = 'trsl_profile_2'
        cProfile.run('tr.solve_internal(0, '+str(duration)+', deltat='+str(deltat)+')', filename)
        import pstats
        p = pstats.Stats(filename)
        p.strip_dirs().sort_stats('cumulative').print_stats()

        tr.dump_results(description)