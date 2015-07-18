'''
Created on 19.02.2015

@author: MJS
'''
import os.path
import cPickle as pkl
import csv
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess  # @UnresolvedImport

datadir = "../data"

plotkinfile = r'Plotkin_calc_IP_output.txt'       
plotkinfilepath = os.path.join(datadir, plotkinfile)
plotkin_ips = {}

with open(plotkinfilepath, mode='r') as infile:
    reader = csv.reader(infile)
    reader.next() # skip header line
    for rows in reader:
        try:
            plotkin_ips[rows[0]] = float(rows[5].replace(',','.'))
        except:
            print "found missing value, skipping..."
plotkin_ips = {key: plotkin_ips[key] for key in plotkin_ips if plotkin_ips[key]>1e-60}
#print plotkin_ips
pkl.dump(plotkin_ips, open("../parameters/init_rates_plotkin.p", "wb"))

# http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002866
# 10.1371/journal.pcbi.1002866
# PLoS Comput Biol 9(1): e1002866. 
stansfield_ips = pkl.load(open("../parameters/init_rates_stansfield.p"))
#print stansfield_ips

plotkeys = stansfield_ips.viewkeys() & plotkin_ips.viewkeys() # intersection
print len(plotkeys)

xs = [stansfield_ips[key] for key in plotkeys]
ys = [plotkin_ips[key] for key in plotkeys] 

fig =  plt.figure()                

ax = fig.add_subplot(111)
ax.grid(True, linestyle = '-', color = '0.75')
#ax.set_xlim([1, 10000])
#ax.set_ylim([1, 10000])
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Stansfield initiation probabilities')
ax.set_ylabel('Plotkin initiation probabilities')

scat = plt.scatter(xs, ys)
scat.set_alpha(0.5)

points = zip(xs, ys)
sorted_points = sorted(points)
new_xs = [point[0] for point in sorted_points]
new_ys = [point[1] for point in sorted_points]
ylowess = lowess(new_ys, new_xs)[:,1]
plt.plot(new_xs, ylowess, 'orange', linewidth=4)

plt.show()

