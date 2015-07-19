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
import random as ran

import numpy.random as npr
import cPickle as pkl

import TRSL
import MRNA
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
    'gug': 'cac', 'gcg': 'ugc', 'gag': 'cuc', 'ggg': 'ucc'
}

''' source: http://nar.oxfordjournals.org/content/suppl/2011/04/23/gkr300.DC1/Supplemental_File_S2.pdf
            http://nar.oxfordjournals.org/content/39/15/6705
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
    def __init__(self, mRNAs, gene_library, decay_constants=None, nribo=200000, proteome=col.Counter({})):
        super(TRSL_spec, self).__init__(nribo, proteome)
        self._tRNA = col.Counter({i: tRNA_types[i]['abundancy'] for i in tRNA_types})
        self.mRNAs = mRNAs
        self.modeldict = {}
        self.decay_constants = decay_constants
        self.initialize_modeldict(mRNAs, gene_library)

    def initialize_modeldict(self, mRNAs, gene_library):
        '''
        Create the modeldict required for the WCM (TODO:)
        The tRNAs, mRNAs and proteins are added with additional indices for separate identification in the loops below
        '''
        self.modeldict['name'] = 'TRSL:_discrete_translation_gene_specific'
        self.modeldict['vars'] = ['ribosomes', 'ATP', 'AMP', 'GTP', 'GDP']
        self.modeldict['initvars'] = col.Counter({'GTP': self.GTP, 'GDP': self.GDP, 'ATP': self.ATP, 'AMP': self.AMP,
                                                  'ribosomes': self.ribo_free})
        self.modeldict['pars'] = []  # TODO: populate pars
        self.modeldict['sp_annotations'] = {'ribosomes': 'GO:0005840',
                                            'GTP': 'CHEBI:15996', 'GDP': 'CHEBI:17552', 'ATP': 'CHEBI:15422',
                                            'AMP': 'CHEBI:16027'}
        self.modeldict['sp_compartment'] = {'ribosomes': 'cytosol', 'GTP': 'cytosol', 'GDP': 'cytosol',
                                            'ATP': 'cytosol', 'AMP': 'cytosol'}
        self.modeldict['units'] = {'ribosomes': 1, 'GTP': 1, 'GDP': 1, 'ATP': 1, 'AMP': 1}

        # add tRNAs:
        for tRNA_num in tRNA_types:
            # the distinction between bound and free tRNA exists only at module level
            self.modeldict['vars'].extend(['tRNA_' + str('%02d' % (tRNA_num,))])
            self.modeldict['initvars']['tRNA_' + str('%02d' % (tRNA_num,))] = tRNA_types[tRNA_num]['abundancy']
            self.modeldict['sp_annotations']['tRNA_' + str('%02d' % (tRNA_num,))] = 'CHEBI:17843_' + str(
                '%02d' % (tRNA_num,))
            # ersetzen durch CHEBIs aus http://www.ebi.ac.uk/ebisearch/search.ebi?db=smallMolecules&t=tRNA
            # oder http://www.ebi.ac.uk/chebi/searchId.do;010D9AC7FDDC72158F86B943C40AD04A?chebiId=CHEBI:2651
            self.modeldict['sp_compartment']['tRNA_' + str('%02d' % (tRNA_num,))] = 'cytosol'
            self.modeldict['units']['tRNA_' + str('%02d' % (tRNA_num,))] = "dimensionless"

        # add mRNAs:
        for mRNA in mRNAs:
            # every mRNA molecule is a separate species because they might have difference ribosomal occupancies
            self.modeldict['vars'].extend(['mRNA_' + str('%02d' % (mRNA.index,))])
            self.modeldict['initvars'][
                'mRNA_' + str('%02d' % (mRNA.index,))] += 1  # this is possible for Counters, not for dicts
            self.modeldict['sp_annotations']['mRNA_' + str('%02d' % (mRNA.index,))] = 'CHEBI:33699_' + str(
                '%02d' % (mRNA.index,))
            self.modeldict['sp_compartment']['mRNA_' + str('%02d' % (mRNA.index,))] = 'cytosol'
            self.modeldict['units']['mRNA_' + str('%02d' % (tRNA_num,))] = 1
            # print mRNA, mRNA.index, mRNA.sequence

        # add proteins:
        for gene in gene_library:
            self.modeldict['vars'].extend(
                ['protein_' + str(gene)])  # FIXME: the proteins should NOT be different, hence labelled by the sequence
            self.modeldict['initvars']['protein_' + str(gene)] = 0
            self.modeldict['sp_annotations']['protein_' + str(gene)] = 'CHEBI:36080_' + str(gene)
            self.modeldict['sp_compartment']['protein_' + str(gene)] = 'cytosol'
            self.modeldict['units']['protein_' + str(gene)] = 1

        self.modeldict['com_annotations'] = {'cytosol': 'GO:0005829'}
        self.modeldict['solver'] = self.solve_internal
        self.modeldict['timerange'] = self.timerange
        self.modeldict['timecourses'] = self.timecourses

    def diffuse_ribosomes_to_initiation_site(self, mRNA, deltat):
        """Perform Poisson experiment to test how many ribosomes make it initiation site and try to attach."""
        if self.ribo_free > 0:
            # k = npr.binomial(self.ribo_free, mRNA.init_rate*deltat, 1)[0]  # number of ribosomes that diffuse to the initiation site during deltat
            k = npr.poisson(
                self.ribo_free * mRNA.init_rate * deltat)  # number of ribosomes that diffuse to the initiation site during deltat
            # log.debug('update_initiation: %s ribosomes out of %s diffused to init site at mRNA %s with probability %s',
            # k, self.ribo_free, mRNA.index, self.init_rate*deltat)
            for i in range(k):  # currently k>1 will not attach k ribosomes, TODO:
                if not mRNA.first_position_occupied():
                    # log.debug("update_initiation: found mRNA with free first position")
                    if self.GTP > 0 and self.ATP > 0:
                        if mRNA.attach_ribosome_at_start():
                            # log.debug("update_initiation: attaching ribosome at start of mRNA %s", mRNA.index)
                            self.ribo_bound += 1
                            self.ribo_free -= 1
                            self.GTP -= 1
                            self.GDP += 1
                            self.ATP -= 2
                            self.AMP += 2

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

    def update_initiation(self, deltat, mRNA):
        # log.info('update_initiation: starting')
        # log.debug('update_initiation: found mRNA %s', mRNA)
        self.diffuse_ribosomes_to_initiation_site(mRNA, deltat)

    def update_elongation(self, deltat, mRNA):
        # log.info("update_elongation: starting")
        # while a change occurs:
        # update all empty ribosomes by tRNA diffusion
        # if possible:
        #   all occupied ribosomes move by one step
        #   after the move they lose bound tRNA
        # halve time interval and continue
        change_flag = True
        available_time = deltat
        while change_flag:  # while there is a change in tRNA or ribosome position
            change_flag = self.fill_empty_ribosomes(mRNA, available_time)  # if a tRNA bound this becomes True
            change_flag = change_flag or self.elongate_mRNA(mRNA)  # if a ribosome translocated this becomes True
            available_time *= 0.5

    def fill_empty_ribosomes(self, mRNA, deltat):
        """Walk through every empty ribosome and try to diffuse the required tRNA into the site."""
        change_occurred = False
        empty_ribos = [key for key in mRNA.ribosomes if mRNA.ribosomes[key] is None]  # TODO: test if termination position must be excluded here
        for ribo_pos in empty_ribos:
            thiscodon = mRNA.sequence[ribo_pos: ribo_pos + 3]
            if thiscodon in stopcodons:
                continue
            print ribo_pos
            print "thiscodon:", thiscodon
            print "codon_anticodon[thiscodon]:", codon_anticodon[thiscodon]
            #import time; time.sleep(50.0 / 1000.0)
            required_tRNA_type = anticodon_index[codon_anticodon[thiscodon]]  # index of anticodon corresponding to next codon in mRNA
            tRNA_diffusion_probability = self.elong_rate * deltat * wobble[thiscodon]
            failure_probability = (1 - tRNA_diffusion_probability) ** self.tRNA_free[required_tRNA_type]
            randomnumber = ran.random()
            # can also try Poisson approximation if faster
            # log.debug("update_initiation: failure probability is %s at mRNA position 0", failure_probability)
            success = not (randomnumber < failure_probability)  # this means the required tRNA type has diffused to the site
            if success:
                # log.debug('update_initiation: matching tRNA diffused to initiation site')
                if not self.insert_tRNA(mRNA, ribo_pos, required_tRNA_type):
                    log.warning("elongate_mRNA: unsuccessful attempt to insert tRNA")
                else:
                    # log.debug("elongate_mRNA: successful attempt to insert tRNA")
                    change_occurred = True
        return change_occurred

    def elongate_mRNA(self, mRNA):
        # log.debug("update_elongation: ribosomes on this mRNA are: %s", mRNA.ribosomes)
        change_occurred = False
        occupied_ribos = [key for key in mRNA.ribosomes if mRNA.ribosomes[key]]
        for ribo_pos in occupied_ribos:
            present_pos = ribo_pos
            available_nucleotides = max(mRNA.find_max_free_range(present_pos) - 3 * MRNA.cr, 0)
            # all empty ribosomes may get occupied by a tRNA
            thiscodon = mRNA.sequence[present_pos: present_pos + 3]
            if thiscodon in stopcodons:
                log.warning("elongate_mRNA: encountered stop codon during elongation step")
                continue
            else:
                change_occurred = self.elongate_one_step(mRNA, present_pos)
        return change_occurred

    def elongate_one_step(self, mRNA, current_pos):
        '''
        attempts to elongate the protein on mRNA by one AA at current_pos
        stops if steric hindrance by another ribosome , or end of mRNA is encountered
        '''
        free_range = mRNA.find_max_free_range(current_pos)
        # log.debug("elongate_one_step: found free range of %s nts downstream of %s", free_range, current_pos)
        codons = min(1, free_range / 3)  # integer division on purpose
        codons = min((mRNA.length - current_pos) / 3, codons)  # cannot translate behind end of mRNA
        # log.debug("elongate_one_step: free range of %s nts, trying to elongate by %s codons", free_range, codons)
        if self.GTP >= codons and codons == 1:
            # log.debug("elongate_one_step: possible to translocate by %s codons", codons)
            # elongation: release tRNA
            previous_type = mRNA.ribosomes[current_pos]  # type to be released at ribo_pos
            # log.debug("elongate_one_step: mRNA.ribosomes = %s", mRNA.ribosomes)
            # log.debug("elongate_one_step: self.tRNA_bound = %s", self.tRNA_bound)
            # translocation: move ribosome
            mRNA.translocate_ribosome(current_pos, by=3)
            self.release_tRNA(mRNA, current_pos+3, previous_type)
            # translocation: elongate proteinlength
            self.protein_length += 1
            self.GTP -= 1
            self.GDP += 1
            return True
        else:
            log.warning("elongate_one_step: not possible: not enough GTP or no codon")
            return False
            # log.debug("elongate_one_step: ribosomes: tRNA is now %s", mRNA.ribosomes)
            # log.debug("elongate_one_step: protein length is now %s", self.proteinlength)

    def update_protein_decay(self, deltat):
        if self.decay_constants:
            # log.info("update_protein_decay: starting")
            for gene in self.proteins:
                # print gene, self.decay_constants[gene]
                self.proteins[gene] *= 1 - deltat * self.decay_constants[
                    gene]  # percentage of surviving proteins (non-integer)
        else:
            # log.warning("update_protein_decay: skipping protein decay")
            pass

    def update_processes(self, deltat):
        for mRNA in self.mRNAs:
            self.update_termination(mRNA)
            self.update_initiation(deltat, mRNA)
            self.update_elongation(deltat, mRNA)
            self.update_protein_decay(deltat)


if __name__ == "__main__":
    log.basicConfig(level=log.DEBUG, format='%(message)s', stream=sys.stdout)

    examplesequence_1 = "ggg uuu uca uca uuu gag gac gau gua ggg uuu uca uca uuu gag gac gau gua ggg uuu uca uca uuu gag gac gau gua uaa".replace(
        ' ', '')
    examplesequence_2 = "aug aaa cug ccc gag ggg uuu uca uca uuu gag gac aaa cug ccc gag ggg uuu uca uca uuu gag gac aaa cug ccc gag ggg uuu uca uca uaa".replace(
        ' ', '')
    '''
    # demo configuration 1: 1 mRNA
    gene_library = {1: examplesequence_1}
    mRNAs = [MRNA_specific.mRNA_spec(index=0, sequence=examplesequence_1, geneID=1)]

    # demo configuration 2: 3 mRNAs of 2 genes
    # all genes that exist and their indices
    gene_library = {
                        1: examplesequence_1,
                        2: examplesequence_2
                        }

    # how many mRNAs there are of a given type
    mRNA_abundancies = {
                        1: 1,
                        2: 2
                        }
    '''

    # demo configuration 3: (part of full-scale) simulation
    # load sequences, transcriptome and initiation rates from pickle file
    mRNA_abundancies = pkl.load(open("../parameters/transcriptome.p", "rb"))
    gene_library = pkl.load(open("../parameters/orf_coding.p", "rb"))
    init_rates = pkl.load(open("../parameters/init_rates_plotkin.p", "rb"))

    mRNAs = []
    counter = 0
    for gene in gene_library.keys()[0:100]:  # shortened for demo purposes
        print "found gene:", gene
        if gene in mRNA_abundancies and gene in init_rates:
            print "abundancies and initiation probabilities available for gene:", gene
            for instance in range(mRNA_abundancies[gene]):
                mRNAs.append(MRNA_specific.mRNA_spec(index=counter, sequence=gene_library[gene], geneID=gene,
                                                     init_rate=init_rates[gene]))  # do not just multiply the list
                counter += 1
    print "built gene library, next: run TRSL_spec."

    tr = TRSL_spec(mRNAs, gene_library, nribo=200000)

    duration = 80.0
    tr.solve_internal(0.0, duration, deltat=1.0)
    # tr.inspect()
