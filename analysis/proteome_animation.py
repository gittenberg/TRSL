'''
Created on 10.01.2015

@author: MJS
'''
import os.path
import cPickle as pkl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from Bio import SeqIO

# get proteome_model and transcriptome data from pickle files 
proteome_model = pkl.load(open("results_20150204_1731_0600s.p", "rb"))
print proteome_model['description']
transcriptome = pkl.load(open("transcriptome.p", "rb"))

# get ORF data data from pickle files 
datadir = "../data"
sequencefile = "orf_coding.fasta"
sequencefile = os.path.join(datadir, sequencefile)
handle = open(sequencefile, "rU")
record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
handle.close()
orf_genomic_dict = {key:str(record_dict[key].seq).lower().replace('t', 'u') for key in record_dict}

plotkeys = (proteome_model["timecourses"].viewkeys() & transcriptome.viewkeys()) & orf_genomic_dict.viewkeys() # intersection; after 300 s there were 681 plotkeys (i.e. finished proteins)
print len(plotkeys)

xs = [transcriptome[key] for key in plotkeys]
ys = [proteome_model["timecourses"][key][0] for key in plotkeys]
nframes = int(proteome_model["duration"])

def _update_plot(i, fig, scat):
    ys = [proteome_model["timecourses"][key][i] for key in plotkeys]
    scat.set_offsets(([[xs[n], ys[n]] for n in range(len(xs))]))
    #print('Frames: %d' %i)
    return scat,

fig =  plt.figure()                

x = xs
y = ys[0]
size = [(len(orf_genomic_dict[key])**1.8)/2000. for key in plotkeys]

ax = fig.add_subplot(111)
ax.grid(True, linestyle = '-', color = '0.75')
ax.set_xlim([1, 10000])
ax.set_ylim([1, 10000])
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('transcript count')
ax.set_ylabel('protein count')

#scat = plt.scatter(xs, ys, c = xs)
scat = plt.scatter(xs, ys, s=size)
scat.set_alpha(0.5)

anim = animation.FuncAnimation(fig, _update_plot, fargs = (fig, scat),
                               frames = nframes, interval=25, repeat_delay=3000, blit=True)

# Set up formatting for the movie files

Writer = animation.writers['ffmpeg']
mywriter = Writer(fps=10, metadata=dict(artist='MJS'), bitrate=1800)
anim.save('translation_flat_initiation.mp4', writer=mywriter)

plt.show()

