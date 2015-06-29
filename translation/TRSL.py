"""
protein translation as a probabilistic, discrete model (specific version)

This module generates an dynamical system of polysome translating proteins.
The proteins are generic (not sequence-specific).

Required steps and energy consumption:

1. initiation (cost: 1 GTP) # TODO: is the probability/rate in the paper including the tRNA or just the ribosome??
1.1 ribosome attachment (first subunit)
1.2 AA-tRNA binding at P site. We assume the AA is activated elsewhere, so no ATP cost
1.3 ribosome attachment (second subunit)
   (mRNA activation:    1 ATP -> ADP, 
   initiation complex:  1 ATP -> ADP, 1 GTP -> GDP)
   We currently merge 1.1, 1.2, 1.3 into one binomial experiment

2. elongation
2.1 AA-tRNA binding at A site (cost: 1*length GTP -> GDP)
2.2 peptide bond formation and peptide elongation (no cost for bond)
2.3 translocation P>E, A>P (cost: 1*length GTP -> GDP)
2.4 release tRNA from E
    We carry out 2.4 (release) before 2.1 (binding) so we do not have to distinguish E, P, A sites

3. termination (cost: 1 GTP -> GDP)

Total energy balance per protein:
2 ATP -> 2 ADP
(2 + 2*length) GTP -> (2 + 2*length) GDP 
"""

import sys
import cProfile
import pstats
import numpy as np 
import numpy.random as npr
import random as ran
import collections as col
import math
import logging as log

import MRNA

class TRSL:
    '''
    class representing a translational network
    '''
    def __init__(self, nribo=200000, proteome=col.Counter({})):
        '''
        initializes the parameters of the translation process
        '''
        log.info("__init__: initializing TRSL")

        # Parameters
        ##################################################################################################################################
        self.types_tRNA = 42               # number of types of tRNAs
        V = 4.2e-17                        # m^3 # cell volume  
        n_genes = 3795                     # number of protein-coding genes in the experiment
        n_tRNA = 3300000                   # number of tRNAs # http://nar.oxfordjournals.org/content/suppl/2011/04/23/gkr300.DC1/Supplemental_File_S2.pdf gives 3000000
             
        avogadro = 6.022e23
        lambda_tRNA = 1.5e-8               # m # characteristic length tRNA
        D_tRNA = 8.42e-11                  # m^2/s # diffusion coefficient of tRNA
        tau_tRNA = lambda_tRNA**2/6/D_tRNA # s # char. time for tRNA
        num_pos_tRNA = V/lambda_tRNA**3    # number of discrete positions for tRNA
        lambda_ribo = 3e-8                 # m # characteristic length tRNA
        D_ribo = 3e-13                     # m^2/s # diffusion coefficient of ribosomes
        tau_ribo = lambda_ribo**2/6/D_ribo # s # char. time for ribosomes
        num_pos_ribo = V/lambda_ribo**3    # number of discrete positions for ribosomes
        p_init = math.sqrt(3.5e-6*0.115)   # initiation probability at mRNA 5' end # we choose the geometric mean btw the lowest and highest possible value
        competition = 7.78e-4              # tRNA competition coefficient

        # Initial values
        ##################################################################################################################################
        self.n_mRNA = 6                    # 60000 # number of mRNAs 
        self.ribo_free = nribo             # 200000; number of ribosomes not bound to mRNA # http://bionumbers.hms.harvard.edu/search.aspx?task=searchbytrmorg&log=y&unt=y&trm=ribosomes/cell&org=%
             
        self.GTP = 1e2*avogadro*V          # GTP molecules (made up)
        self.GDP = 0                       # GDP molecules
        self.ATP = 1e2*avogadro*V          # ATP molecules (made up)
        self.AMP = 0                       # AMP molecules
        #self.AA = 1e2*avogadro*V          # amino acids (made up)

        self.timerange = []
        self.timecourses = {}
        self.ribo_bound = 0                # number of ribosomes bound to mRNA
        self.tRNA = col.Counter({i:int(0.5+n_tRNA/self.types_tRNA) for i in range(1, self.types_tRNA+1)}) # TODO: experimental values are in Ingolia (2009)
        self.tRNA_free = col.Counter({i:int(self.tRNA[i]) for i in range(1, self.types_tRNA+1)}) # tRNA not bound to ribosomes
        self.tRNA_bound = self.tRNA - self.tRNA_free                                             # tRNA bound to ribosomes

        self.mRNAs = [MRNA.MRNA(index=gene) for gene in [ran.randint(1, n_genes) for k in range(self.n_mRNA)]] # randomized gene expressions
        self.proteins = proteome  # contains protein IDs and counts not including polypeptides in statu nascendi 
        self.proteinlength = sum(self.proteins.values())
        
        self.init_rate = p_init/tau_ribo/num_pos_ribo               # 8.2e-07 s^-1 (yeast)
        self.elong_rate = competition/tau_tRNA/num_pos_tRNA         # 0.0001 s^-1  (yeast)

        self.modeldict = {'name': "TRSL:_discrete_translation",
                          'vars': ["protein", "ribos._bound", "ribos._free", "tRNA_bound", "tRNA_free", "ATP", "AMP", "GTP", "GDP"],
                          'initvars': {"protein": 0, "ribos._bound": 0, "ribos._free": self.ribo_free, "tRNA_bound": 0, "tRNA_free": sum(self.tRNA_free.values()), 'GTP': self.GTP, 'GDP': self.GDP, 'ATP': self.ATP, 'AMP': self.AMP},
                          'pars': [],
                          'sp_annotations': {"protein": "CHEBI:36080",    # generic protein
                                             #"amino_acid": "CHEBI:15705", # generic amino acid
                                             "ribos._bound": "GO:0042788", 
                                             "ribos._free":  "GO:0005840", 
                                             "tRNA_bound": "CHEBI:17843_b", 
                                             "tRNA_free": "GO:0005564", # generic tRNA; http://www.ebi.ac.uk/chebi/searchId.do;010D9AC7FDDC72158F86B943C40AD04A?chebiId=CHEBI:2651 lists some others 
                                             'GTP': 'CHEBI:15996', 'GDP': 'CHEBI:17552', 'ATP': 'CHEBI:15422', 'AMP': 'CHEBI:16027'},
                          'sp_compartment': {"protein": 'cytosol', "ribos._bound": 'cytosol', "ribos._free": 'cytosol', "tRNA_bound": 'cytosol', "tRNA_free": 'cytosol', 'GTP': 'cytosol', 'GDP': 'cytosol', 'ATP': 'cytosol', 'AMP': 'cytosol'},
                          'units': {"protein": 1, "ribos._bound": 1, "ribos._free": 1, "tRNA_bound": 1, "tRNA_free": 1, 'GTP': 1, 'GDP': 1, 'ATP': 1, 'AMP': 1},
                          'com_annotations': {'cytosol': 'GO:0005829'},
                          'solver': self.solve_internal,
                          'timerange': self.timerange,
                          'timecourses': self.timecourses
                          } # TODO: work in progress

    def __getitem__(self, i):
        """
        This allows to address the modeldict like trsl['vars']
        """
        return self.modeldict[i]
            
    def inspect(self):
        """
        Print all dictionaries in TRSL
        """
        print "----------------------------------------------------------"
        for key in sorted(self.__dict__):
            if key != "modeldict":
                print key, ":"
                print self.__dict__[key]
                print "----------------------------------------------------------"
            else:
                print key, ":"
                for subkey in self.__dict__[key]:
                    print "\t\t", subkey, ":"
                    #print "\t\t", self.__dict__[key][subkey]
                    print "----------------------------------------------------------"
    
    def dump_results(self, description='results'):
        results = {}
        results['proteome'] = self.proteins
        results['transcriptome'] = self.mRNAs
        results['timerange'] = self.timerange
        results['timecourses'] = self.timecourses
        results["description"] = description
        import time; now = time.strftime("%Y%m%d_%H%M", time.gmtime())
        results["time_stamp"] = now
        results["n_ribosomes"] = self.ribo_bound + self.ribo_free
        results["n_tRNA"] = sum(self.tRNA.values())
        duration = self.timerange[-1] - self.timerange[0]
        results["duration"] = duration
        from cPickle import dump; dump(results, open("../results/results_"+now+"_"+str(int(duration)).zfill(4)+"s.p", "wb"))
        print description

    def insert_tRNA(self, mRNA, pos, tRNA_type):
        '''
        inserts a tRNA of type tRNA_type at position pos
        returns True iff successful
        '''
        if pos in mRNA.ribosomes and self.tRNA_free[tRNA_type]>=1 and not mRNA.ribosomes[pos]:
            # there has to be a ribosome at pos and there has to be tRNA of that type available and there cannot be a tRNA yet at pos on the position
            #log.debug("insert_tRNA: inserting tRNA %s on mRNA %s at position %s", tRNA_type, self.mRNAs.index(mRNA), pos)
            self.tRNA_free[tRNA_type] -= 1
            self.tRNA_bound[tRNA_type] += 1
            mRNA.ribosomes[pos] = tRNA_type  # tRNA now bound
            #log.debug("insert_tRNA: ribosomes: tRNA on mRNA %s are now %s", self.mRNAs.index(mRNA), mRNA.ribosomes)
            success = True
        elif not pos in mRNA.ribosomes:
            log.warning("insert_tRNA: cannot insert tRNA at %s: no ribosome present", pos)
            success = False
        elif mRNA.ribosomes[pos]!=None:
            log.warning("insert_tRNA: cannot insert tRNA at %s: tRNA already bound", pos)
            success = False
        elif self.tRNA_free[tRNA_type]<1:
            #log.warning("insert_tRNA: cannot insert tRNA type %s at %s: not enough free tRNA", tRNA_type, pos)
            #log.warning("insert_tRNA: tRNA_type is: %s", tRNA_type)
            #log.warning("insert_tRNA: bound tRNA is: %s", self.tRNA_bound)
            success = False
        else:
            log.warning("insert_tRNA: failed for unspecified reason")
            success = False
        return success

    def release_tRNA(self, mRNA, pos, tRNA_type):
        '''
        releases a tRNA molecule from mRNA
        '''
        if mRNA.ribosomes[pos]==tRNA_type and self.tRNA_bound[tRNA_type]>0:
            self.tRNA_bound[tRNA_type] -= 1
            self.tRNA_free[tRNA_type] += 1
            mRNA.ribosomes[pos] = None
            #log.debug("release_tRNA: successfully released tRNA %s from pos %s", tRNA_type, pos)
            success = True
        else:
            log.warning("release_tRNA: failed to release tRNA from pos %s", pos)
            log.warning("release_tRNA: mRNA.ribosomes[pos] == %s", mRNA.ribosomes[pos])
            log.warning("release_tRNA: self.tRNA_bound[tRNA_type] == %s", self.tRNA_bound[tRNA_type])
            success = False
        return success
    
    def update_initiation(self, deltat):
        # log.info("update_initiation: starting")
        for mRNA in self.mRNAs:
            #log.debug("update_initiation: found mRNA %s", j)
            if self.ribo_free>0:
                k = npr.binomial(self.ribo_free, self.init_rate*deltat, 1)[0]  # number of ribosomes that diffuse to the initiation site during deltat
                #log.debug("update_initiation: %s ribosomes diffused to init site at mRNA %s", k, mRNA.index)
                for i in range(k):  # currently k>1 will not attach k ribosomes, perhaps after parallelization?
                    if not mRNA.first_position_occupied():
                        #log.debug("update_initiation: found mRNA with free first position")
                        if self.GTP>0 and self.ATP>0:
                            if mRNA.attach_ribosome_at_start():
                                #log.debug("update_initiation: attaching ribosome at start of mRNA %s", mRNA.index)
                                self.ribo_bound += 1
                                self.ribo_free -= 1
                                self.GTP -= 1
                                self.GDP += 1
                                self.ATP -= 2
                                self.AMP += 2
                                init_type = ran.choice(self.tRNA.keys())  # type to be inserted at pos==0
                                if not self.insert_tRNA(mRNA, 0, init_type):
                                    log.warning("update_initiation: unsuccessful attempt to insert tRNA")
                            else:
                                log.warning("update_initiation: unsuccessful attempt to attach ribosome")
                        else:
                            log.warning("update_initiation: no GTP or no ATP")
                    else: 
                        #log.warning("update_initiation: unsuccessful attempt to attach ribosome: first position occupied")
                        pass
            else: 
                #log.warning("update_initiation: no free ribosomes left")
                pass
            
    def elongate_while_possible(self, mRNA, k, current_pos):
        '''
        attempts to elongate the protein on mRNA by at most k AAs starting at current_pos
        stops if steric hindrance by another ribosome , or end of mRNA is encountered
        '''
        free_range = mRNA.find_max_free_range(current_pos)
        #log.debug("elongate_while_possible: found free range of %s nts downstream of %s", free_range, current_pos)
        codons = min(k, free_range/3) # number k of available tRNAs and sterically free codons limit elongation # integer division
        codons = max(0, codons) # not negative
        codons = min((mRNA.length-current_pos)/3, codons) # cannot translate behind end of mRNA
        #log.debug("elongate_while_possible: %s tRNAs, free range of %s nts, trying to elongate by %s codons", k, free_range, codons)
        if self.GTP >= codons:
            if codons>0:
                #log.debug("elongate_while_possible: possible to translocate by %s codons", codons)
                # elongation: release tRNA
                previous_type = mRNA.ribosomes[current_pos] # type to be released at ribo_pos, could be made sequence dependent
                self.release_tRNA(mRNA, current_pos, previous_type)
                # translocation: move ribosome
                mRNA.translocate_ribosome(current_pos, by=3*codons)
                # bind AA-tRNA
                #last_type = ran.choice(range(1, self.types_tRNA+1)) # type to be inserted at ribo_pos
                last_type = ran.choice(self.tRNA.keys()) # type to be inserted at current_pos
                #log.debug("elongate_while_possible: last position was %s, attempting to insert tRNA at position %s", current_pos, current_pos+3*codons)
                self.insert_tRNA(mRNA, current_pos+3*codons, last_type) # try to insert AA-tRNA in the ribosome
                # translocation: elongate proteinlength
                self.proteinlength += codons
                self.GTP -= codons
                self.GDP += codons
        else:
            log.warning("elongate_while_possible: not enough GTP")
        #log.debug("elongate_while_possible: ribosomes: tRNA is now %s", mRNA.ribosomes)
        #log.debug("elongate_while_possible: protein length is now %s", self.proteinlength)

    def update_elongation(self, deltat):
        for mRNA in self.mRNAs:
            #for ribo_pos in npr.permutation(mRNA.ribosomes.keys()):  # random permutation of the ribosomes # slow!
            #log.debug("update_elongation: ribosomes on this mRNA are: %s", mRNA.ribosomes)
            for ribo_pos in mRNA.ribosomes.keys():
                present_tRNA_type = mRNA.ribosomes[ribo_pos]
                k = npr.binomial(self.tRNA_free[present_tRNA_type], self.elong_rate*deltat, 1)[0]  # number of tRNA_free[present_tRNA_type] that diffuse to the ribosome during deltat
                log.debug("update_elongation: %s tRNAs out of %s of type %s diffused to elongation site %s at mRNA %s", k, self.tRNA_free[present_tRNA_type], present_tRNA_type, ribo_pos, mRNA.index)
                current_pos = ribo_pos
                self.elongate_while_possible(mRNA, k, current_pos)

    def update_termination(self):
        log.info("update_termination: starting")
        for mRNA in self.mRNAs:
            if self.GTP >= 1:
                if mRNA.termination_condition():
                    #log.debug("update_termination: mRNA.ribosomes = %s", mRNA.ribosomes)
                    release_pos = max(mRNA.ribosomes.keys())
                    release_type = mRNA.ribosomes[release_pos]
                    #log.debug("update_termination: terminating translation at mRNA %s, release position %s, release type %s", mRNA.index, release_pos, release_type)
                    self.release_tRNA(mRNA, release_pos, release_type)
                    mRNA.detach_ribosome(release_pos)
                    self.ribo_bound -= 1
                    self.ribo_free += 1
                    if not mRNA.geneID in self.proteins:
                        self.proteins[mRNA.geneID] = 1 # add first protein of type mRNA.geneID
                    else:
                        self.proteins[mRNA.geneID] += 1 # add another protein of type mRNA.geneID
                    self.GTP -= 1
                    self.GDP += 1
            else:
                log.warning("update_termination: not enough GTP")
    
    def update_processes(self, deltat):
        self.update_termination()
        self.update_initiation(deltat)
        self.update_elongation(deltat)
        #self.update_protein_decay(deltat)
        
    def solve_internal(self, start, end, deltat):
        '''
        solves TRSL for the interval [start, end[, iterating through several steps
        '''
        log.info("solve: simulation from %s to %s", start, end)

        fieldnames = ["protein", "ribos._bound", "ribos._free", "tRNA_bound", "tRNA_free", "ATP", "AMP", "GTP", "GDP"]
        self.timecourses = {fieldname: [] for fieldname in fieldnames}

        for tRNA_type in self.tRNA_free:
            self.timecourses["tRNA_free_"+str(tRNA_type).zfill(2)] = []
        
        self.timerange = np.arange(start, end, deltat)
        for time in self.timerange:
            log.info("########################################################################################################")
            log.info("solve: time: %s", time)
            log.info("########################################################################################################")

            self.update_processes(deltat)
            log.info("solve_internal: self.proteins = %s", self.proteins)
            log.info("solve_internal: protein length:  %s", self.proteinlength)
            log.info("solve_internal: bound ribosomes: %s", self.ribo_bound)
            log.info("solve_internal: free ribosomes:  %s", self.ribo_free)
            log.info("solve_internal: bound tRNA:      %s", sum(self.tRNA_bound.values()))
            log.info("solve_internal: free tRNA:       %s", sum(self.tRNA_free.values()))
            fieldvalues = [self.proteinlength, self.ribo_bound, self.ribo_free, sum(self.tRNA_bound.values()), sum(self.tRNA_free.values()), self.ATP, self.AMP, self.GTP, self.GDP]

            # update everything except proteins and specific tRNA_free
            for fieldname, fieldvalue in zip(fieldnames, fieldvalues):
                self.timecourses[fieldname].append(fieldvalue)

            # now update proteins
            for gene_id in self.proteins:
                if gene_id in self.timecourses: # if there is already protein of this ID 
                    # then append the protein count
                    self.timecourses[gene_id].append(self.proteins[gene_id])
                else:                           # if this is the first time protein of this ID terminated
                    # create zeros for the first (time-start)/deltat values
                    self.timecourses[gene_id] = [0] * int((time - start) / deltat)
                    # only then append the protein count
                    self.timecourses[gene_id].append(self.proteins[gene_id])

            # now update specific tRNA_free
            for tRNA_type in self.tRNA_free:
                self.timecourses["tRNA_free_"+str(tRNA_type).zfill(2)].append(self.tRNA_free[tRNA_type])
                #log.info("solve_internal: tRNA_free type %s: %s molecules", tRNA_type, self.tRNA_free[tRNA_type])

        #from time import gmtime, strftime; now = strftime("%Y%m%d_%H%M%S", gmtime())
        #import os.path; pickle.dump({'trange': self.timerange, "timecourses": self.timecourses}, open(os.path.join("..", now+"_timecourses_TRSL.pkl"), "wb"))

        return self.proteins, self.mRNAs


if __name__=="__main__":
    log.basicConfig(level=log.INFO, format='%(message)s', stream=sys.stdout)
    trsl=TRSL()
    '''
    trsl.solve_internal(0.0, 60.0, deltat=1.0)
    '''
    # Profiling:
    cProfile.run('trsl.solve_internal(0.0, 60.0, deltat=1.0)', 'trsl_profile')
    p=pstats.Stats('trsl_profile')
    p.strip_dirs().sort_stats('cumulative').print_stats('TRSL')
