'''
Created on 27.04.2014

Makes plots based on previous run of the TRSL model

@author: MJS
'''
import pickle as pkl
import os
import matplotlib.pyplot as plt
import numpy as np
import math
#import os.path

def extract_ribos_mrnas(f):
    fsplit = f.split("_")
    ribo_index = fsplit.index("ribos") - 1 # where to find the number of ribosomes
    mRNA_index = fsplit.index("mRNAs") - 1 # where to find the number of mRNAs
    ribos = int(fsplit[ribo_index])
    mRNAs = int(fsplit[mRNA_index])
    return ribos, mRNAs

def get_data(f):
    data = pkl.load(file(f))
    file.close(file(f))
    return data

def plotTimecourses(trange, timecourses, log_y_scale=True, autooff=False,
                              titlestring=None):
    """
    Plot the time course of the species in a certain module.

    Parameters
    ----------
    trange : list
        time points for the x axis
    timecourses : dict
        dict of all available species data points
    one_plot : Boolean
        defaul True, should all plots be in one graph
    log_y_scale : Boolean
        default True, should the y axis use a logarithm scale
    """
    fig, ax = plt.subplots()
    if log_y_scale:
        plt.yscale('log')
    if titlestring:
        ax.set_title(titlestring)
    else:
        ax.set_title('Click on legend line to toggle line on/off')
    line = {}
    # we will set up a dict mapping legend line to orig line, and enable picking on the legend line
    lined = {}
    maxlength = max([len(timecourses[speciesname]) for speciesname in timecourses])
    for speciesname in timecourses:
        if len(timecourses[speciesname])==maxlength:
            line[(speciesname)] = ax.plot(trange, timecourses[speciesname], lw=1, label=speciesname)

    leg = ax.legend(loc='best', fancybox=True, shadow=False, prop={'size': 10})
    leg.get_frame().set_alpha(0.4)
    
    for legline, origline in zip(leg.get_lines(), ax.lines):
        legline.set_picker(5)  # 5 pts tolerance
        lined[legline] = origline
        if autooff:
            legline.set_alpha(0.2)
            origline.set_visible(False)
    
    def onpick(event):
        # on the pick event, find the orig line corresponding to the
        # legend proxy line, and toggle the visibility
        legline = event.artist
        origline = lined[legline]
        vis = not origline.get_visible()
        origline.set_visible(vis)
        # Change the alpha on the line in the legend so we can see what lines have been toggled
        if vis:
            legline.set_alpha(1.0)
        else:
            legline.set_alpha(0.2)
        fig.canvas.draw()
    
    fig.canvas.mpl_connect('pick_event', onpick)
    plt.savefig(titlestring+'.png')
    #plt.show()

def plot3D(x, y, z, title, logscale):
    from mpl_toolkits.mplot3d import Axes3D  # @UnresolvedImport
    
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
    ax.w_xaxis.set_ticklabels(row_names)
    ax.w_yaxis.set_ticklabels(column_names)
    ax.set_xlabel('mRNAs')
    ax.set_ylabel('ribosomes')
    ax.set_zlabel(title)
    
    plt.show()

if __name__=="__main__":
    pkl_dir = os.listdir(os.path.join(os.getcwd(), "../results"))
    print pkl_dir
    x, y = [], []
    peptide_bonds, proteinrate = [], []

    for f in pkl_dir:
        if f.endswith(".p"):
            #print f
            ribos, mRNAs = extract_ribos_mrnas(f)
            #print ribos, mRNAs 
            data = get_data(os.path.join(os.getcwd(), "analysis_data", f))
            trange = data["trange"]
            timecourses = data["timecourses"]
            #print timecourses.keys()
            
            proteinfinal = timecourses["protein"][-1]
            proteininital = timecourses["protein"][0]
            x.append(ribos)
            y.append(mRNAs)
            peptide_bonds.append(proteinfinal)
            diffquot = (proteinfinal - proteininital)/(trange[-1] - trange[0])
            proteinrate.append(diffquot)
            
            ribobound = timecourses["ribos._bound"][-1] # latest value
            if ribos == 200000:
                print 1/(ribobound/float(mRNAs)/1250*3) # codon distance of ribosomes

            displayed_timecourses = {key: timecourses[key] for key in ["protein"]}
            titlestring = str(ribos)+" ribosomes, "+str(mRNAs)+" mRNAs"
            plotTimecourses(trange, displayed_timecourses, log_y_scale=False, autooff=False, titlestring=titlestring)
            
    x = [6, 60, 600, 6000]
    y = [2, 20, 200, 2000]
    print peptide_bonds
    plot3D(x, y, peptide_bonds, title="log(peptide bonds + 1)", logscale=True)
    print proteinrate
    plot3D(x, y, proteinrate, title="log(protein rate + 1)", logscale=True)