'''
Created on 18.04.2015
Show how to connect to keypress events

key 'j' moves one step forward
key 'k' moves one step backward

@author: MJS
'''
import sys
import matplotlib.pyplot as plt
import cPickle as pkl

prot = pkl.load(open("results_20150603_1414_7200s_Plotkin.p", "rb"))

print prot['description']

timerange = prot['timerange']
timecourses = prot['timecourses']

keys = timecourses.keys()
n = len(keys)
i = 0
fig, ax = plt.subplots()
this_protein_ID = keys[i]
xs, ys = timerange, timecourses[this_protein_ID]
line1, = ax.plot(xs, ys, 'g-')
plt.ylim([0, max(ys)])
ax.set_xlabel('time [s]')
ax.set_ylabel(this_protein_ID)

def press(event):
    global i
    sys.stdout.flush()
    direction = 0
    if event.key == 'j':
        direction = 1
    elif event.key == 'k':
        direction = -1
    next_protein_ID = keys[(i + direction) % n]
    xs, ys = timerange, timecourses[next_protein_ID]
    line1.set_xdata(xs)
    line1.set_ydata(ys)
    plt.ylim([0, max(ys)])
    ax.set_ylabel(next_protein_ID)
    fig.canvas.draw()
    i += direction

fig.canvas.mpl_connect('key_press_event', press)
plt.show()
