import sys
import cPickle as pkl
import logging as log

from translation import MRNA_specific, TRSL_specific

log.basicConfig(level=log.DEBUG, format='%(message)s', stream=sys.stdout)

# transcriptome from Weinberg (private email Shah, 2015)
transcriptome = pkl.load(open("../parameters/transcriptome_shah.p", "rb"))

# initiation rates from Weinberg (private email Shah, 2015)
init_rates = pkl.load(open("../parameters/init_rates_plotkin.p", "rb"))

# exome
exome = pkl.load(open("../parameters/orf_coding.p", "rb"))

# initiation rates scaling
scaling_factors = [0.01, 0.1, 1.0, 10.0, 100.0]

for scale in scaling_factors:
    scaled_init_rates = {key: scale * init_rates[key] for key in init_rates}

    # description
    description = 'updated Shah transcriptome, full exome, no decay, updated initiation rates according to Shah, initiation rates scaled by ' + str(scale)
    print description

    genes = list(set(exome) & set(transcriptome) & set(scaled_init_rates))
    print "found %s genes in common." % len(genes)

    mRNAs = []
    counter = 0
    for gene in genes:
        for transcript in range(transcriptome[gene]):
            # create mRNA with scaled initiation rate
            mRNAs.append(MRNA_specific.mRNA_spec(index=counter, sequence=exome[gene], geneID=gene, ribosomes={},
                                                 init_rate=scaled_init_rates[gene]))  # do not just multiply the list
            counter += 1
    print "built gene library, next: run TRSL_spec."

    tr = TRSL_specific.TRSL_spec(mRNAs, exome, decay_constants=None, nribo=200000, detail=True)

    print "found {} mRNAs".format(len(mRNAs))
    print "found {} genes".format(len(exome))
    print "found {} tRNA molecules".format(sum(tr._tRNA.values()))
    print "solving..."

    duration = 1800.0
    tr.solve_internal(0.0, duration, deltat=0.2)
