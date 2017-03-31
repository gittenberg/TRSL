"""
protein translation as a probabilistic, discrete model (specific version)

This module generates an dynamical system of polysome translating proteins.
The proteins are sequence-specific.

Module inherits from generic module TRSL.

@author: martin
"""

import sys
import logging as log
import collections as col
import numpy as np

import numpy.random as npr
import cPickle as pkl

import TRSL
import MRNA_specific


#############################################################################################################
# auxiliary data
#############################################################################################################

'''
these codons are on the mRNA from 5' to 3'
'''
genetic_code = {
    'uuu': 'F', 'ucu': 'S', 'uau': 'Y', 'ugu': 'C',
    'uuc': 'F', 'ucc': 'S', 'uac': 'Y', 'ugc': 'C',
    'uua': 'L', 'uca': 'S', 'uaa': '*', 'uga': '*',  # '*'==stop
    'uug': 'L', 'ucg': 'S', 'uag': '*', 'ugg': 'W',
    'cuu': 'L', 'ccu': 'P', 'cau': 'H', 'cgu': 'R',
    'cuc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
    'cua': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
    'cug': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
    'auu': 'I', 'acu': 'T', 'aau': 'N', 'agu': 'S',
    'auc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
    'aua': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
    'aug': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
    'guu': 'V', 'gcu': 'A', 'gau': 'D', 'ggu': 'G',
    'guc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
    'gua': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
    'gug': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'
}

'''
source: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC147333/
we use 5'-3' convention for codon and anticodon;
we ignore special nucleotides i and psi and use a and u instead;
some simplification by making this a 1:1 relationship
'''

codon_anticodon = {
    'uuu': 'gaa', 'ucu': 'aga', 'uau': 'gua', 'ugu': 'gca',
    'uuc': 'gaa', 'ucc': 'aga', 'uac': 'gua', 'ugc': 'gca',
    'uua': 'uaa', 'uca': 'uga', 'uaa': '*', 'uga': '*',  # '*'==stop
    'uug': 'caa', 'ucg': 'cga', 'uag': '*', 'ugg': 'cca',
    'cuu': 'gag', 'ccu': 'agg', 'cau': 'gug', 'cgu': 'acg',
    'cuc': 'gag', 'ccc': 'agg', 'cac': 'gug', 'cgc': 'acg',
    'cua': 'uag', 'cca': 'ugg', 'caa': 'uug', 'cga': 'acg',
    'cug': 'uag', 'ccg': 'ugg', 'cag': 'cug', 'cgg': 'ccg',
    'auu': 'aau', 'acu': 'agu', 'aau': 'guu', 'agu': 'gcu',
    'auc': 'aau', 'acc': 'agu', 'aac': 'guu', 'agc': 'gcu',
    'aua': 'uau', 'aca': 'ugu', 'aaa': 'uuu', 'aga': 'ucu',
    'aug': 'cau', 'acg': 'cgu', 'aag': 'cuu', 'agg': 'ccu',
    'guu': 'aac', 'gcu': 'agc', 'gau': 'guc', 'ggu': 'gcc',
    'guc': 'aac', 'gcc': 'agc', 'gac': 'guc', 'ggc': 'gcc',
    'gua': 'uac', 'gca': 'ugc', 'gaa': 'uuc', 'gga': 'ucc',
    'gug': 'cac', 'gcg': 'ugc', 'gag': 'cuc', 'ggg': 'ccc'
}

''' source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3159466/bin/supp_39_15_6705__index.html
            https://academic.oup.com/nar/article/39/15/6705/1022014/The-role-of-tRNA-and-ribosome-competition-in
    anticodons are in 5'-3' direction (same convention as above)'''
tRNA_types = {
    1: {'anticodon': 'ugc', 'abundancy': 55351},  # reverse complement the anticodon to look it up
    2: {'anticodon': 'agc', 'abundancy': 121771},
    3: {'anticodon': 'ucu', 'abundancy': 121771},
    4: {'anticodon': 'ccu', 'abundancy': 11070},
    5: {'anticodon': 'ccg', 'abundancy': 11070},
    6: {'anticodon': 'acg', 'abundancy': 66421},
    7: {'anticodon': 'guu', 'abundancy': 110701},
    8: {'anticodon': 'guc', 'abundancy': 177122},
    9: {'anticodon': 'gca', 'abundancy': 44280},
    10: {'anticodon': 'uug', 'abundancy': 88561},
    11: {'anticodon': 'cug', 'abundancy': 11070},
    12: {'anticodon': 'uuc', 'abundancy': 154982},
    13: {'anticodon': 'ucc', 'abundancy': 33210},
    14: {'anticodon': 'ccc', 'abundancy': 22140},
    15: {'anticodon': 'gcc', 'abundancy': 177122},
    16: {'anticodon': 'gug', 'abundancy': 77491},
    17: {'anticodon': 'uau', 'abundancy': 22140},
    18: {'anticodon': 'aau', 'abundancy': 143911},
    19: {'anticodon': 'uag', 'abundancy': 33210},
    20: {'anticodon': 'gag', 'abundancy': 11070},
    21: {'anticodon': 'uaa', 'abundancy': 77491},
    22: {'anticodon': 'caa', 'abundancy': 110701},
    23: {'anticodon': 'uuu', 'abundancy': 77491},
    24: {'anticodon': 'cuu', 'abundancy': 154982},
    25: {'anticodon': 'cau', 'abundancy': 55351},  # 26 does not seem to exist
    27: {'anticodon': 'gaa', 'abundancy': 110701},
    28: {'anticodon': 'agg', 'abundancy': 22140},
    29: {'anticodon': 'ugg', 'abundancy': 110701},
    30: {'anticodon': 'gcu', 'abundancy': 33210},
    31: {'anticodon': 'uga', 'abundancy': 33210},
    32: {'anticodon': 'aga', 'abundancy': 121771},
    33: {'anticodon': 'cga', 'abundancy': 11070},
    34: {'anticodon': 'ugu', 'abundancy': 44280},
    35: {'anticodon': 'agu', 'abundancy': 121771},
    36: {'anticodon': 'cgu', 'abundancy': 11070},
    37: {'anticodon': 'cca', 'abundancy': 66421},
    38: {'anticodon': 'gua', 'abundancy': 88561},
    39: {'anticodon': 'uac', 'abundancy': 22140},
    40: {'anticodon': 'aac', 'abundancy': 154982},
    41: {'anticodon': 'cac', 'abundancy': 22140},
    42: {'anticodon': '*', 'abundancy': 18000},  # termination factor
    43: {'anticodon': 'cuc', 'abundancy': 22140}
}

stopcodons = ['uaa', 'uga', 'uag']
anticodon_index = {tRNA_types[i]['anticodon']: i for i in tRNA_types}


#############################################################################################################
# auxiliary functions
#############################################################################################################

def tRNA_type_at_position(sequence, pos):
    return anticodon_index[codon_anticodon[sequence[pos:pos + 3]]]


def chunker(seq, size):
    '''
    generator that takes a sequence and returns substrings of given size
    '''
    return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))


def translate_mRNA(sequence):
    '''
    Converts an mRNA sequence to a correspondingly translated protein sequence
    '''
    codons = chunker(sequence, 3)
    protein = ""
    for nxt in codons:
        if nxt not in stopcodons:
            protein += genetic_code[nxt]
        else:
            return protein
    log.warning("translate_mRNA: WARNING: no stop codon found, end of mRNA reached. Returning protein")
    return protein


def complement(s):
    ''' complementary RNA'''
    basecomplement = {'a': 'u', 'c': 'g', 'g': 'c', 'u': 'a'}
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)


def revcom(s):
    ''' reverse complementary RNA'''
    if s == "*":
        return False
    return complement(s[::-1])


wobble = {codon: 1.000 if codon_anticodon[codon] is not '*' and codon == revcom(codon_anticodon[codon]) else 0.625 for
          codon in codon_anticodon}


#############################################################################################################
# TRSL class definition
#############################################################################################################

class TRSL_spec(TRSL.TRSL):

    def __init__(self, mRNAs, gene_library, decay_constants=None, nribo=200000, proteome=col.Counter({}), detail=False):
        super(TRSL_spec, self).__init__(nribo, proteome, detail)
        self._tRNA = col.Counter({i: tRNA_types[i]['abundancy'] for i in tRNA_types})
        self.mRNAs = mRNAs
        self.n_mRNA = len(self.mRNAs)
        self.modeldict = {}
        self.decay_constants = decay_constants


    def diffuse_ribosomes_to_initiation_site(self, mRNA, deltat, time):
        """Perform Poisson experiment to test how many ribosomes make it initiation site and try to attach."""
        if self.ribo_free > 0:
            # k = npr.binomial(self.ribo_free, mRNA.init_rate*deltat, 1)[0]  # number of ribosomes that diffuse to the initiation site during deltat
            k = npr.poisson(
                self.ribo_free * mRNA.init_rate * deltat)  # number of ribosomes that diffuse to the initiation site during deltat
            # log.debug('update_initiation: %s ribosomes out of %s diffused to init site at mRNA %s with probability %s', k, self.ribo_free, mRNA.index, self.init_rate*deltat)
            for i in range(k):  # currently k>1 will not attach k ribosomes, TODO:
                if not mRNA.ribosomes or not mRNA.first_position_occupied():
                    # log.debug("update_initiation: found mRNA with free first position")
                    if self.GTP > 0 and self.ATP > 0:
                        if mRNA.attach_ribosome_at_start():
                            # log.debug("update_initiation: attaching ribosome at start of mRNA %s", mRNA.index)
                            self.ribo_bound += 1
                            self.ribo_free -= 1
                            self.GTP -= 1
                            self.GDP += 1
                            self.ATP -= 1
                            self.AMP += 1
                            if not mRNA.tic:  # no time measurement ongoing on this mRNA
                                mRNA.tic = time
                                mRNA.toc = len(mRNA.ribosomes) + 1  # number of ribos + 1 to countdown to end of time measurement
                        else:
                            log.warning("update_initiation: unsuccessful attempt to attach ribosome")
                    else:
                        log.warning("update_initiation: no GTP or no ATP")
                else:
                    # log.warning("update_initiation: unsuccessful attempt to attach ribosome: first position occupied")
                    pass
        else:
            # log.warning("update_initiation: no free ribosomes left")
            pass

    def update_initiation(self, deltat, mRNA, time):
        # log.info('update_initiation: starting')
        # log.debug('update_initiation: found mRNA %s', mRNA)
        self.diffuse_ribosomes_to_initiation_site(mRNA, deltat, time)  # tic = True if an initiation occurred

    def fill_empty_ribosomes(self, mRNA, deltat):
        """Walk through every empty ribosome and try to diffuse the required tRNA into the site."""
        change_occurred = False
        # empty ribosomes on this particular mRNA
        # the ribosome at position 0 is filled during the initiation process (not modelled)
        empty_ribos = [key for key in mRNA.ribosomes if mRNA.ribosomes[key] is None and key!=0]
        # TODO: test if termination position must be excluded here
        for ribo_pos in empty_ribos:
            thiscodon = mRNA.sequence[ribo_pos: ribo_pos + 3]
            if thiscodon in stopcodons:
                continue

            required_tRNA_type = anticodon_index[codon_anticodon[thiscodon]]  # index of anticodon corresponding to next codon in mRNA
            tRNA_diffusion_probability = self.elong_rate * deltat * wobble[thiscodon]
            failure_probability = (1 - tRNA_diffusion_probability) ** self.tRNA_free[required_tRNA_type]
            randomnumber = np.random.ranf()  # TODO: try Poisson approximation if faster
            # log.debug("fill_empty_ribosomes: failure probability is %s at mRNA position 0", failure_probability)
            success = not (randomnumber < failure_probability)  # this means the required tRNA type has diffused to the site
            if success:
                # log.debug('fill_empty_ribosomes: matching tRNA diffused to initiation site')
                if not self.insert_tRNA(mRNA, ribo_pos, required_tRNA_type):
                    log.warning("elongate_mRNA: unsuccessful attempt to insert tRNA")
                else:
                    # log.debug("elongate_mRNA: successful attempt to insert tRNA")
                    change_occurred = True
        return change_occurred

    def elongate_mRNA(self, mRNA):
        """translocates all ribosomes on mRNA by one step"""
        # log.debug("elongate_mRNA: ribosomes on this mRNA are: %s", mRNA.ribosomes)
        occupied_ribos = [key for key in mRNA.ribosomes if mRNA.ribosomes[key] is not None or key==0] # the first codon is always occupied by tRNA^Met_i
        for ribo_pos in occupied_ribos:
            present_pos = ribo_pos
            thiscodon = mRNA.sequence[present_pos: present_pos + 3]  # TODO: is this redundant because elongate_one_step or translocate_ribosome are testing the same?
            if thiscodon in stopcodons:
                #log.warning("elongate_mRNA: encountered stop codon during elongation step, not elongating this one")
                continue
            else:
                self.elongate_one_step(mRNA, present_pos)

    def update_protein_decay(self, deltat):
        if self.decay_constants:
            # log.info("update_protein_decay: starting")
            for gene in self.proteins:
                # print gene, self.decay_constants[gene]
                self.proteins[gene] *= 1 - deltat * self.decay_constants[gene]  # percentage of surviving proteins (non-integer)
        else:
            # log.warning("update_protein_decay: skipping protein decay")
            pass

    def update_processes(self, deltat, time):
        for mRNA in self.mRNAs:
            if mRNA.ribosomes != {}:  # mRNAs without ribosomes do not need the following
                self.update_elongation(deltat, mRNA)
                self.update_termination(mRNA, time)
            self.update_initiation(deltat, mRNA, time)
            #self.update_protein_decay(deltat)


if __name__ == "__main__":
    log.basicConfig(level=log.DEBUG, format='%(message)s', stream=sys.stdout)

    conf = {
           'exome': pkl.load(open("../parameters/orf_coding.p", "rb")),
           'transcriptome': pkl.load(open("../parameters/transcriptome_plotkin.p", "rb")),
           'init_rates': pkl.load(open("../parameters/init_rates_plotkin.p", "rb")),
           'decay_constants': pkl.load(open("../parameters/decay_constants.p", "rb")),
           'description': 'full transcriptome and exome, with decay, specific best estimate initiation rates according to Plotkin'
           }

    genes = list(set(conf['exome']) & set(conf['transcriptome']) & set(conf['init_rates']) & set(conf['decay_constants']))
    print "found %s genes in common." % len(genes)

    mRNAs = []
    counter = 0
    for gene in genes:
        # print "abundancies and initiation rates available for gene:", gene
        for instance in range(conf['transcriptome'][gene]):
            mRNAs.append(MRNA_specific.mRNA_spec(index=counter, sequence=conf['exome'][gene], geneID=gene, ribosomes={}, init_rate=conf['init_rates'][gene]))  # do not just multiply the list
            counter += 1
    print "built gene library, next: run TRSL_spec."

    description = conf['description']
    print description

    duration = 20.0

    tr = TRSL_spec(mRNAs, conf['exome'], conf['decay_constants'], nribo=20000)

    tr._tRNA = col.Counter({i: tRNA_types[i]['abundancy'] for i in tRNA_types})
    tr._tRNA_free = col.Counter({i: int(tr._tRNA[i]) for i in tRNA_types})  # tRNA not bound to ribosomes
    tr._tRNA_bound = tr._tRNA - tr._tRNA_free  # tRNA bound to ribosomes
    tr.solve_internal(0.0, duration, deltat=1.0)
    tr.dump_results(description='results')

    '''
    # Profiling:
    import cProfile
    cProfile.run('tr.solve_internal(0.0, '+str(duration)+', deltat=1.0)', 'trsl_profile')
    import pstats
    p=pstats.Stats('trsl_profile')
    p.strip_dirs().sort_stats('cumulative').print_stats()
    '''