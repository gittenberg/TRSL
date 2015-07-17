'''
Created on 04.11.2014

@author: martin
'''
import matplotlib.pyplot as plt
import numpy

mrnas = [6, 60, 600, 6000, 60000]
times = [0.082, 0.824, 7.664, 83.700, 449.09]

fig, ax = plt.subplots()
ind = numpy.arange(len(mrnas))
width = 0.35
plt.bar(ind, times, log=True)
ax.set_title('Time required to solve timecourses [s]', fontsize=11)
ax.set_xlabel('number of mRNAs')
ax.set_xticks(ind + width)
ax.set_xticklabels(mrnas)
plt.show()

# test
