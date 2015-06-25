'''
Created on 22.02.2015

@author: martin
'''
import collections as col
import cPickle as pkl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess

print "plotting proteomes..."

# get proteome_exp
'''
import os.path
datadir = "../data"
# http://pax-db.org/#!downloads
#proteomefile = r"4932-Spectral_counting_S.cerevisiae_GPM_Oct_2012.txt"
#proteomefile = r'4932-Spectral_counting_S.cerevisiae_PeptideAtlas_May_2009.txt'
#proteomefile = r'4932-S.cerevisiae_whole_organism-integrated_dataset.txt'       # combined data
proteomefile = r'4932-S.cerevisiae_whole_organism-integrated_dataset.txt'       # combined data
proteomefilepath = os.path.join(datadir, proteomefile)
with open(proteomefilepath, mode='r') as infile:
    lines = infile.readlines()[9:]
    lines = [(line.split()[1], line.split()[2]) for line in lines]
    proteome_exp = dict((line[0].split('.')[1], float(line[1])) for line in lines)
'''
# we use the same proteome as Liebermeister to be consistent:
proteome_exp = pkl.load(open("proteome_Nagaraj.p", "rb"))

# get proteome_model and transcriptome data from pickle files 
proteome_model = pkl.load(open("results_20150204_1552_0600s.p", "rb"))
print proteome_model['description']                          # features of the underlying transcriptome

annotations = pkl.load(open("annotations_Liebermeister.p", "rb"))
functional_categories = col.Counter(annotations.values())
'''
for i, j in enumerate(functional_categories.most_common()):
    print i, j[0]
'''
colors = cm.rainbow(np.linspace(0, 1, len(functional_categories)))  # @UndefinedVariable
colordict = {j[0]: colors[i] for i, j in enumerate(functional_categories.most_common())}
print colordict

plotkeys = proteome_model["timecourses"].viewkeys() & proteome_exp.viewkeys() & annotations.viewkeys() # intersection
print len(plotkeys)

xs = [proteome_exp[key] for key in plotkeys]
ys = [proteome_model["timecourses"][key][-1] for key in plotkeys] # use the last value as a proxy for the proteome_model

fig =  plt.figure()                

ax = fig.add_subplot(111)
ax.grid(True, linestyle = '-', color = '0.75')
#ax.set_xlim([1, 10000])
#ax.set_ylim([1, 10000])
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Protein abundance (measured)')
ax.set_ylabel('Protein abundance (modelled)')


colors = [colordict[annotations[key]] for key in plotkeys]
scat = plt.scatter(xs, ys, c=colors)
scat.set_alpha(0.5)
#plt.legend()

points = zip(xs, ys)
sorted_points = sorted(points)
new_xs = [point[0] for point in sorted_points]
new_ys = [point[1] for point in sorted_points]
ylowess = lowess(new_ys, new_xs)[:,1]
plt.plot(new_xs, ylowess, 'orange', linewidth=4)
'''
#plt.show()

print "plotting proteome vs. transcriptome..."

transcriptome = pkl.load(open("transcriptome.p", "rb"))
plotkeys = proteome_exp.viewkeys() & transcriptome.viewkeys() # intersection; after 300 s there were 681 plotkeys (i.e. finished proteins)

xs = [proteome_exp[key] for key in plotkeys]
ys = [transcriptome[key] for key in plotkeys]

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
'''
print "done."

plt.show()

