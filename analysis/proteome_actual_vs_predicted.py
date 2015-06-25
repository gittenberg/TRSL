'''
Created on 21.01.2015

@author: martin
'''
import os.path
import cPickle as pkl
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess  # @UnresolvedImport
from scipy.stats.stats import pearsonr   

print "plotting proteomes..."

# get proteome_exp
datadir = "../data"
# http://pax-db.org/#!downloads
#proteomefile = r"4932-Spectral_counting_S.cerevisiae_GPM_Oct_2012.txt"
#proteomefile = r'4932-Spectral_counting_S.cerevisiae_PeptideAtlas_May_2009.txt'
proteomefile = r'4932-S.cerevisiae_whole_organism-integrated_dataset.txt'       # combined data
proteomefilepath = os.path.join(datadir, proteomefile)

with open(proteomefilepath, mode='r') as infile:
    lines = infile.readlines()[9:]
    lines = [(line.split()[1], line.split()[2]) for line in lines]
    proteome_exp = dict((line[0].split('.')[1], float(line[1])) for line in lines)

# get proteome_model
proteome_model = pkl.load(open("results_20150603_1414_7200s_Plotkin.p", "rb"))
print proteome_model['description'] # description of the underlying simulation

# get transcriptome
transcriptome = pkl.load(open("transcriptome.p", "rb"))

plotkeys = proteome_model["timecourses"].viewkeys() & proteome_exp.viewkeys() # intersection; after 300 s there were 681 plotkeys (i.e. finished proteins)
#plotkeys = proteome_model["timecourses"].viewkeys() & proteome_exp.viewkeys() & transcriptome.viewkeys() 
# to put both comparisons (actual vs. predicted and actual vs. transcriptome on the same basis

print len(plotkeys)

xs = [proteome_exp[key] for key in plotkeys]
ys = [proteome_model["timecourses"][key][-1] for key in plotkeys] # use the last value as a proxy for the proteome_model

print "correlation coefficient: {0:.2f}%".format(pearsonr(xs, ys)[0]*100) # [1] is the p-value which is not interesting

fig =  plt.figure()                

ax = fig.add_subplot(111)
ax.grid(True, linestyle = '-', color = '0.75')
ax.set_xlim([1, 10000])
ax.set_ylim([1, 1000])
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Protein abundance (measured)')
ax.set_ylabel('Protein abundance (modelled)')

scat = plt.scatter(xs, ys)
scat.set_alpha(0.5)

points = zip(xs, ys)
sorted_points = sorted(points)
new_xs = [point[0] for point in sorted_points]
new_ys = [point[1] for point in sorted_points]
ylowess = lowess(new_ys, new_xs)[:,1]
plt.plot(new_xs, ylowess, 'orange', linewidth=4)

#plt.show()

# transcriptome

print
print "plotting proteome vs. transcriptome..."

plotkeys = proteome_exp.viewkeys() & transcriptome.viewkeys() # intersection; after 300 s there were 681 plotkeys (i.e. finished proteins)
#plotkeys = proteome_model["timecourses"].viewkeys() & proteome_exp.viewkeys() & transcriptome.viewkeys() 
# to put both comparisons (actual vs. predicted and actual vs. transcriptome on the same basis

xs = [proteome_exp[key] for key in plotkeys]
ys = [transcriptome[key] for key in plotkeys]

print "correlation coefficient: {0:.2f}%".format(pearsonr(xs, ys)[0]*100) # [1] is the p-value which is not interesting

fig2 =  plt.figure()                

ax = fig2.add_subplot(111)
ax.grid(True, linestyle = '-', color = '0.75')
ax.set_xlim([1, 10000])
ax.set_ylim([1, 10000])
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Protein abundance (measured)')
ax.set_ylabel('Transcript abundance (measured)')

scat = plt.scatter(xs, ys)
scat.set_alpha(0.5)

points = zip(xs, ys)
sorted_points = sorted(points)
new_xs = [point[0] for point in sorted_points]
new_ys = [point[1] for point in sorted_points]
ylowess = lowess(new_ys, new_xs)[:,1]
plt.plot(new_xs, ylowess, 'orange', linewidth=4)

print "done."

plt.show()
