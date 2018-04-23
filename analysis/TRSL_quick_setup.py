'''
Created on 23.04.2018

To test the maximum of ribosome efficiency in a quicker setting (fewer transcripts, tRNAs)

@author: MJS
'''

import logging as log
import sys
import collections as col
import datetime
import os
import os.path

from translation import MRNA_specific, TRSL_specific

print "TRSL_vary_ribosomes: starting"
log.basicConfig(level=log.DEBUG, format='%(message)s', stream=sys.stdout)

examplesequence_1 = "ggg uuu uca uca uuu gag gac gau gua aaa ggg uaa".replace(' ', '')
examplesequence_2 = "aug aaa cug ccc gag ggg uuu uca uca uuu gag gac aaa cug ccc gag ggg uuu uca uca uuu gag gac aaa cug ccc gag ggg uuu uca uca uaa".replace(' ', '')

# configuration
exome = {1: examplesequence_1, 2: examplesequence_2}
transcriptome = {1: 200, 2: 300}
init_rates = {1: 0.1, 2: 0.001}
description_short = 'test configuration with 2 genes in 50 transcripts'
tRNA = col.Counter({i: 100 for i in TRSL_specific.tRNA_types})

print "found {} transcripts in transcriptome...".format(sum(transcriptome.values()))

duration = 200.0  # 1200. should be sufficent to saturate
ribonumbers = [1, 2, 3, 5, 8, 10, 20, 30, 50, 80, 100, 200]

# loop over different ribosome counts
for ribonumber in ribonumbers:
    deltat = 0.05
    description = '{}, {} ribosomes, deltat={}s'.format(description_short, ribonumber, deltat)

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
    directory = os.path.join('..', 'results_'+today)
    if not os.path.exists(directory):
        os.makedirs(directory)
    tr.dump_results(description, dirname=directory)
