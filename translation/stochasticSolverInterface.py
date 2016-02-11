

from abc import ABCMeta, abstractmethod

class StochasticSolverInterface(object):
        __metaclass__ = ABCMeta

        def __init__(self):
            pass

        @abstractmethod
        def execute(self, trange, x0 = None):
            """
            :param trange: list of time steps, depending on step size.
                           eg: starttime = 0
                               endtime = 1
                               step size = 0.1
                               --> trange list = [0, 0.1, 0.2, 0.3,....,1.0]
                    x0:    initial state vector as a dictionary. It's not given at first call
                           units is number of particles.
                           eg. x0 = {'ATP': 100, 'NADH':5000, 'GLC':700}

            :return  return stateVector, info
                            stateVector: same type as x0 with last particle numbers
                            info:        status message for solverstep
                                         it's a  dictionary like info = {'message': 'everythings fine'}


            :raise NotImplementedError:
            """
            raise NotImplementedError("function execute(observer)")
