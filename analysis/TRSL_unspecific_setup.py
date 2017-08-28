"""
Produce results for demo charts of unspecific TRSL

Further analysis in: http://localhost:8888/notebooks/workbooks/charts/TRSL_plot_unspecific_timecourses.ipynb
"""
import sys
import logging as log

from translation import MRNA, TRSL

log.basicConfig(level=log.INFO, format='%(message)s', stream=sys.stdout)

"""
Low ribosome scenario, vary number of transcripts
"""

"""
factor = [1, 10, 100, 1000, 10000]

for n in factor:
    print "n = {}, simulating...".format(n)
    trsl = TRSL.TRSL(nribo=20*n)  # 20, 200, 2000, 20000, 200000
    # overwrite number of transcripts:
    trsl.n_mRNA = 6 * n           # 6, 60, 600, 6000, 60000

    trsl.solve_internal(0.0, 300.0, deltat=0.2)
    trsl.dump_results(description='TRSL_unspecific_low_ribosomes_results_{}_transcripts'.format(n))
"""

"""
Realistic ribosome and transcript scenario
"""
trsl = TRSL.TRSL(nribo=200000)
trsl.n_mRNA = 60000
trsl.solve_internal(0.0, 100.0, deltat=0.2)
trsl.dump_results(description='TRSL_unspecific_realistic_results_{}_transcripts'.format(trsl.n_mRNA))
