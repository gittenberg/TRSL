import re

__author__ = 'max'


"""
Take a pickle file produced by the pypy version of numpy and reformat it to the cnumpy format.
"""

def clean(path):
    """
    Substitute all references of _numpypy to the correct counterparts in cnumpy.

    :param path: path
    """
    cleaned = ""
    with open(path) as f:
        cleaned = re.sub("c_numpypy", "cnumpy.core", f.read())
    with open(path, 'w') as f:
        f.write(cleaned)
    return 1
