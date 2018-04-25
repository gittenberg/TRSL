'''
Created on 23.04.2018

To test the maximum of ribosome efficiency in a quicker setting (fewer transcripts, tRNAs)

@author: MJS
'''

import logging as log
import sys
import collections as col
import datetime
import os.path
import itertools

from translation import MRNA_specific, TRSL_specific

print "TRSL_vary_ribosomes: starting"
log.basicConfig(level=log.DEBUG, format='%(message)s', stream=sys.stdout)

examplesequence_1 = "ggg uuu uca uca uuu gag gac gau gua aaa ggg uaa".replace(' ', '')
examplesequence_2 = "aug aaa cug ccc gag ggg uuu uca uca uuu gag gac aaa cug ccc gag ggg uuu uca uca uuu gag gac aaa cug ccc gag ggg uuu uca uca uaa".replace(' ', '')

# configuration
exome = {0: examplesequence_1, 1: examplesequence_2}
transcriptomes = [{0: 2, 1: 4}, {0: 20, 1: 40}, {0: 200, 1: 400}, {0: 4, 1: 2}, {0: 40, 1: 20}, {0: 400, 1: 200}]
init_rates_list = [{0: 0.1, 1: 0.001}, {0: 0.01, 1: 0.001}, {0: 0.001, 1: 0.001}, {0: 0.001, 1: 0.1}, {0: 0.001, 1: 0.01}, {0: 0.001, 1: 0.1}]

scenarios = []
# add description to scenarios:
for transcriptome, init_rates in list(itertools.product(transcriptomes, init_rates_list)):
    description_short = 'test configuration with 2 genes, transcripts {}, {}, init rates {}, {}'.format(transcriptome[0], transcriptome[1], init_rates[0], init_rates[1])
    scenarios.append((transcriptome, init_rates, description_short))
    print scenarios

tRNA = col.Counter({i: 100 for i in TRSL_specific.tRNA_types})

print "found {} transcripts in transcriptome...".format(sum(transcriptome.values()))

duration = 200.0  # 1200. should be sufficent to saturate
ribonumbers = [1, 2, 3, 5, 8, 10, 20, 30, 50, 80, 100, 200, 500]

# loop over scenarios:
for index, (transcriptome, init_rates, description_short) in enumerate(scenarios):
    # loop over different ribosome counts
    for ribonumber in ribonumbers:
        deltat = 0.1
        description = '{}, {} ribosomes, deltat={}s'.format(description_short, ribonumber, deltat)
        print description

        counter = 0
        genes = list(set(exome) & set(transcriptome) & set(init_rates))
        print "found %s genes in common." % len(genes)

        # create mRNAs
        mRNAs = []
        for gene in genes:
            for transcript in range(transcriptome[gene]):
                mRNAs.append(MRNA_specific.mRNA_spec(index=counter, sequence=exome[gene], geneID=gene, ribosomes={}, init_rate=init_rates[gene]))  # do not just multiply the list
                counter += 1

        print "created transcriptome: {}.".format(description)

        tr = TRSL_specific.TRSL_spec(mRNAs, exome, decay_constants=None, nribo=ribonumber, proteome=col.Counter({}), detail=False)
        tr.solve_internal(0, duration, deltat=deltat)
        today = datetime.datetime.today().strftime('%Y%m%d')
        dirname = 'results_' + today + '_' + str(index)
        directory = os.path.join('..', dirname)
        if not os.path.exists(directory):
            os.makedirs(directory)
        tr.dump_results(description, dirname=directory)
