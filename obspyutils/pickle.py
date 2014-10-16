# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================

import cPickle

# ----------------------------------------------------------------------
def pickle(filename, data):

    fout = open(filename, "w")
    pickler = cPickle.Pickler(fout, protocol=-1)
    pickler.dump(data)
    fout.close()

    return


# ----------------------------------------------------------------------
def unpickle(filename):

    fin = open(filename, "r")
    unpickler = cPickle.Unpickler(fin)
    data = unpickler.load()
    fin.close()

    return data


# End of file
