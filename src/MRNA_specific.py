"""
mRNA_specific class definition

This module generates an mRNA object representing an mRNA molecule with a specific code sequence.
It inherits from the generic mRNA class.

Parameters:
      index - the ID of this particular mRNA object in the cell
      length - the length in nucleotides
      geneID - the ID of the gene belonging to this particular mRNA object (there might be more than one mRNA sharing the same geneID)
      sequence - the string of nucleotides (in 'a', 'c', 'g', 'u')
      init_rate - the rate (in s^-1) at which initiation start at this mRNA. The default value of the initiation rate  p_init/tau_ribo/num_pos_ribo is calculated in TRSL.py
"""

import logging as log

import MRNA

class mRNA_spec(MRNA.MRNA):
    def __init__(self, index, sequence, geneID, init_rate=8.2e-07): 
        '''
        initializes one mRNA molecule
        '''
        self.sequence = sequence
        self.init_rate = init_rate        # ORF-specific initiation rate # TODO: use this in TRSL_specific.update_initiation
        MRNA.MRNA.__init__(self, index=index, length=len(self.sequence), geneID=geneID)

    def termination_condition(self):
        '''
        returns True iff a ribosome hits a stop codon on the mRNA
        '''
        if self.ribosomes:
            last_pos = max(self.ribosomes.keys())
            stopcodons = ['uga', 'uaa', 'uag']
            current_codon = self.sequence[last_pos:last_pos+3]
            if current_codon in stopcodons: # if ribosome hits stop codon
                #log.debug("termination_condition: found stop codon %s at position %s", current_codon, last_pos)
                return True
            else:
                return False
        else:
            return False
