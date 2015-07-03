__author__ = 'max'

"""
A number of functions to analyze the results of simulations.
"""

import sys
import matplotlib.pyplot as plt
import cPickle as pkl

i = 0

def plot_timecourse(results, subset=None, interactive=True):
    """
    Plot time-courses for results of simulations. Plot the protein amount for each protein
    over the whole simulation time. The method can run interactively by going through proteins
    by pressing k/j for forward/backward movement.

    :param results: set of results from pickle
    :param subset: list of gene names
    :param interactive: run interactive or plot all in one plot
    :return: None
    """
    print results['description']

    timerange = results['timerange']
    timecourses = results['timecourses']

    keys = timecourses.keys()
    n = len(keys)
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

def plot_polysome_histograms(result, subset=None, single=False):
    """
    Plot histogram of amount of ribosomes on mRNA per gene.

    :param genes: list of gene ids
    :param single:
    :return:
    """
    ribos_per_gene = get_polysome_histogram(result)
    # if subset:
    #     n = len(subset)
    #     if n < 20:
    #         sub.x =
    #         plt.subplot(n, )
    if subset:
        for gene in subset:
            plt.figure()
            plt.hist(ribos_per_gene[gene], normed=True, stacked=True)
    else:
        print "Please select a subset to plot."



def get_polysome_histogram(result):
    """
    Compute histogram of amount of ribosomes on mRNA per gene.

    :param result: result dict
    :return: dict of counters
    """

    polysomes_per_gene = {}
    transcripts = result["transcriptome"]
    for mrna in transcripts:
        if mrna.geneID in polysomes_per_gene:
            polysomes_per_gene[mrna.geneID].append(mrna.ribosomes)
        else:
            polysomes_per_gene[mrna.geneID] = [mrna.ribosomes]

    ribosome_count = {gene_id:
                          [len(p) for p in polysomes_per_gene[gene_id]]
                      for gene_id in polysomes_per_gene}
    return ribosome_count




if __name__ == "__main___":
    steady = pkl.load(open("results/steady_state_20150701_2026_0089s.p", "rb"))
    plot_polysome_histograms(steady, subset=("YJL124C","YHL015W"))
    plt.figure()
    run_off = pkl.load(open("results/glucose_starvation_after_steady_20150701_2027_0009s.p", "rb"))
    plot_polysome_histograms(run_off, subset=("YJL124C","YHL015W"))
    plt.show()

    steady_ribos = get_polysome_histogram(steady)
    average_ribos = {x:float(sum(steady_ribos[x]))/len(steady_ribos[x]) for x in steady_ribos}
    max_ribos = {x:max(steady_ribos[x]) for x in steady_ribos}