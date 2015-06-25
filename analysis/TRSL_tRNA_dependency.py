'''

Evaluate the TRSL_specific model by looping over tRNA attenuation and duration

@author: MJS
'''

import TRSL_specific
import MRNA_specific
import sys
import logging as log
import collections as col
import matplotlib.pyplot as plt
import numpy as np
import math

log.basicConfig(level=log.DEBUG, format='%(message)s', stream=sys.stdout)
savedict = {}
savedict['data'] = {}
#alphas = [0.25, 0.5, 0.75, 1.0]
alphas = np.arange(0.0, 1.1, 0.1)

description = "Constant sequence, all tRNAs are variable"
examplesequence = "ggg ggg ggg ggg ggg ggg ggg ggg ggg ggg ggg uaa".replace(' ', '') # uniform sequence
attenuation_factors = {i:alphas for i in TRSL_specific.tRNA_types}                   # all tRNAs are attenuated uniformly

description = "Constant sequence, no tRNA knocked out"
examplesequence = "ggg ggg ggg ggg ggg ggg ggg ggg ggg ggg ggg uaa".replace(' ', '') # uniform sequence
attenuation_factors = {i:[1 for alpha in alphas] for i in TRSL_specific.tRNA_types}  # all tRNAs are not attenuated 

description = "Constant sequence, essential tRNA knocked out"
examplesequence = "ggg ggg ggg ggg ggg ggg ggg ggg ggg ggg ggg uaa".replace(' ', '') # uniform sequence
attenuation_factors = {i:[1 for alpha in alphas] for i in TRSL_specific.tRNA_types}  # all tRNAs are not attenuated... 
attenuation_factors[13] = [0 for alpha in alphas]                                    #... except the relevant one is knocked out

description = "Constant sequence, essential tRNA attenuated"
examplesequence = "ggg ggg ggg ggg ggg ggg ggg ggg ggg ggg ggg uaa".replace(' ', '') # uniform sequence
attenuation_factors = {i:[1 for alpha in alphas] for i in TRSL_specific.tRNA_types}  # all tRNAs are not attenuated... 
attenuation_factors[13] = alphas                                                     #... except the relevant one is attenuated

description = "Mixed sequence, all tRNAs are variable"
examplesequence = "gug ucu ugc uau ucg agu gca ucc aga uca uac aag ugu aua cuu aug cgu gua gac cug cga cuc cua auu gau ggu guc \
                   guu gaa uuu ugg cgc auc cgg gcc ggg gga ggc gag acg uug uua ccg uuc acc aca aaa gcg aac ccu cag gcu agc cau \
                   aau agg acu cac caa cca ccc uaa".replace(' ', '')
attenuation_factors = {i:alphas for i in TRSL_specific.tRNA_types}                   # all tRNAs are attenuated uniformly

description = "Mixed sequence, no tRNA knocked out"
examplesequence = "gug ucu ugc uau ucg agu gca ucc aga uca uac aag ugu aua cuu aug cgu gua gac cug cga cuc cua auu gau ggu guc \
                   guu gaa uuu ugg cgc auc cgg gcc ggg gga ggc gag acg uug uua ccg uuc acc aca aaa gcg aac ccu cag gcu agc cau \
                   aau agg acu cac caa cca ccc uaa".replace(' ', '')
attenuation_factors = {i:[1 for alpha in alphas] for i in TRSL_specific.tRNA_types}  # all tRNAs are not attenuated 

description = "Mixed sequence, essential tRNA (2nd in sequence) knocked out"
examplesequence = "gug ucu ugc uau ucg agu gca ucc aga uca uac aag ugu aua cuu aug cgu gua gac cug cga cuc cua auu gau ggu guc \
                   guu gaa uuu ugg cgc auc cgg gcc ggg gga ggc gag acg uug uua ccg uuc acc aca aaa gcg aac ccu cag gcu agc cau \
                   aau agg acu cac caa cca ccc uaa".replace(' ', '')
attenuation_factors = {i:[1 for alpha in alphas] for i in TRSL_specific.tRNA_types}  # all tRNAs are not attenuated... 
attenuation_factors[32] = [0 for alpha in alphas]                                    #... except the second one required is knocked out

description = "Mixed sequence, essential tRNA (2nd in sequence) attenuated"
examplesequence = "gug ucu ugc uau ucg agu gca ucc aga uca uac aag ugu aua cuu aug cgu gua gac cug cga cuc cua auu gau ggu guc \
                   guu gaa uuu ugg cgc auc cgg gcc ggg gga ggc gag acg uug uua ccg uuc acc aca aaa gcg aac ccu cag gcu agc cau \
                   aau agg acu cac caa cca ccc uaa".replace(' ', '')
attenuation_factors = {i:[1 for alpha in alphas] for i in TRSL_specific.tRNA_types}  # all tRNAs are not attenuated... 
attenuation_factors[32] = alphas                                                     #... except the second one required is knocked out

# change the network's tRNA content:

savedict['experiment'] = description
savedict['sequence'] = examplesequence
durations = [15.0, 30.0, 60.0, 120.0, 240.0, 480.0]
peptide_bonds, proteins = [], []

for duration in durations:
    print "#############################################################################################################"
    print "duration =", duration
    print "#############################################################################################################"

    for k, alpha in enumerate(alphas):
        print "#############################################################################################################"
        print "alpha =", alpha
        print "#############################################################################################################"
        
        # 1 mRNA
        mRNAs = [MRNA_specific.mRNA_spec(index=0, sequence=examplesequence, geneID=1)]
        gene_library = {1: examplesequence}
        
        tr = TRSL_specific.TRSL_spec(mRNAs, gene_library)
        tr.tRNA = col.Counter({i:int(TRSL_specific.tRNA_types[i]['abundancy']*attenuation_factors[i][k]) for i in TRSL_specific.tRNA_types})
        #tr.tRNA = col.Counter({i:TRSL_specific.tRNA_types[i]['abundancy'] for i in TRSL_specific.tRNA_types})
        #tr.tRNA[13] = 0 # we modify a relevant tRNA
        tr.tRNA_free = col.Counter({i:int(TRSL_specific.tRNA_types[i]['abundancy']*attenuation_factors[i][k]) for i in TRSL_specific.tRNA_types})
        #tr.tRNA_free = col.Counter({i:TRSL_specific.tRNA_types[i]['abundancy'] for i in TRSL_specific.tRNA_types})
        #tr.tRNA_free[13] = 0 # we modify a relevant tRNA
        tr.tRNA_bound = col.Counter({i:0 for i in TRSL_specific.tRNA_types})
    
        print "tRNA =", tr.tRNA
        print "tRNA_free =", tr.tRNA_free
        print "tRNA_bound =", tr.tRNA_bound
    
        tr.solve_internal(0.0, duration, deltat=1.0)
        
        peptide_bonds.append(tr.proteinlength)
        proteins.append(tr.proteins)
        savedict['data'][(duration, alpha)] = {"peptide_bonds":tr.proteinlength, "proteins":tr.proteins}

print "#############################################################################################################"
#print "peptide_bonds:", peptide_bonds
#print "proteins:     ", proteins
print savedict

# pickle everything with timestamp
from time import gmtime, strftime
now = strftime("%Y%m%d_%H%M%S", gmtime())
pickle_filename = now + "_mRNAs_timecourses_TRSL.pkl"
import pickle
import os.path
pickle.dump(savedict, open(os.path.join(pickle_filename), "wb"))

def plot3D(x, y, z, xtitle, ytitle, ztitle, header, logscale):
    from mpl_toolkits.mplot3d import Axes3D
    
    row_names = x
    column_names = y
    if logscale:
        z = [math.log(elem+1) for elem in z]
    
    fig = plt.figure()
    ax = Axes3D(fig)
    
    lx = len(x)            # Work out matrix dimensions
    ly = len(y)
    xpos = np.arange(0, lx, 1)    # Set up a mesh of positions
    ypos = np.arange(0, ly, 1)
    xpos, ypos = np.meshgrid(xpos, ypos)
    xpos = xpos.flatten()   # Convert positions to 1D array
    ypos = ypos.flatten()
    zpos = np.zeros(lx*ly)
    
    dx = np.ones_like(zpos)
    dy = dx.copy()
    dz = np.array(z)
    print dz
    
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b')
    
    #sh()
    posxticks = np.arange(len(x))+0.5
    plt.xticks(posxticks, x)
    posyticks = np.arange(len(y))+0.5
    plt.yticks(posyticks, y)
    #ax.w_xaxis.set_ticklabels(range(1, 11))
    ax.w_xaxis.set_ticklabels(row_names)
    ax.w_yaxis.set_ticklabels(column_names)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    ax.set_zlabel(ztitle)
    fig.suptitle(header, fontsize=10)    
    plt.show()

x = alphas
y = durations
plot3D(x, y, peptide_bonds, xtitle='tRNA', ytitle='duration [s]', ztitle='pept. bonds', header=savedict['experiment'], logscale=False)
#plot3D(x, y, peptide_bonds, xtitle='realisation', ytitle='duration [s]', ztitle='pept. bonds', header=savedict['experiment'], logscale=False)
