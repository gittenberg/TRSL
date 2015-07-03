'''
Created on 20.01.2015

@author: martin
'''
import cPickle as pkl
import matplotlib.pyplot as plt

# get proteome_model and transcriptome data from pickle files 
#proteome_model = pkl.load(open("results_20150204_1552_0600s.p", "rb"))
proteome_model = pkl.load(open("../results/glucose_starvation_after_steady_20150701_2027_0009s.p", "rb"))
print proteome_model['description']

timerange = proteome_model['timerange']
fieldnames = ["protein", "ribos._bound", "ribos._free", "tRNA_bound", "tRNA_free", "ATP", "AMP", "GTP", "GDP"]
#for tRNA_type in TRSL_specific.tRNA_types:
#    fieldnames.append("tRNA_free_"+str(tRNA_type).zfill(2))
    
#plt.rcParams.update({'font.size': 15})
for key in fieldnames:
    fig =  plt.figure()                

    print key
    print proteome_model['timecourses'][key]
    plt.plot(timerange, proteome_model['timecourses'][key], 'r-')
    plt.ylim(ymin=0)
    plt.title(key)

fig =  plt.figure()                
ax = fig.add_subplot(111)
ax.set_ylim([0, 220000])
ax.set_xlim([0, 100])
ax.set_xlabel('time [s]')
ax.set_ylabel('bound ribosomes')
plt.plot(timerange, proteome_model['timecourses']['ribos._bound'], 'r-')
plt.plot(timerange, [200000]*len(timerange), 'r--')

plt.show()

'''
fig =  plt.figure()                
ax = fig.add_subplot(111)
#ax.set_ylim([0, 0.05])
ax.set_xlim([0, 20])
import matplotlib.ticker as mtick
fmt = '%.0f%%' # Format you want the ticks, e.g. '40%'
yticks = mtick.FormatStrFormatter(fmt)
ax.yaxis.set_major_formatter(yticks)

ax.set_xlabel('time [s]')
ax.set_ylabel('bound tRNA')
plt.plot(timerange, np.array(proteome_model['timecourses']['tRNA_bound'])/29667.88, 'r-') # total mRNA
plt.ylim(ymin=0)
plt.show()

fig =  plt.figure()                
ax = fig.add_subplot(111)
tRNA_20_final = proteome_model['timecourses']['tRNA_free_20']
proteome_model_unspec = pkl.load(open("results_20150204_1731_0600s.p", "rb"))
tRNA_20_unspec = proteome_model_unspec['timecourses']['tRNA_free_20']
p1 = plt.plot(timerange, tRNA_20_final, 'r-')
p2 = plt.plot(timerange, tRNA_20_unspec, 'b-')
ax.set_ylabel('free tRNA, anticodon GAG')
ax.set_xlabel('time [s]')

plt.legend([p1[0], p2[0]], ['specific initiation', 'unspecific initiation'])

plt.show()

'''
