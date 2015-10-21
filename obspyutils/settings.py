# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================


#-----------------------------------------------------------------------
def merge(a, b):
    """
    Merge dictionaries. Destination is b.
    """
    from copy import deepcopy

    if isinstance(b, dict) and isinstance(a, dict):
        a_and_b = a.viewkeys() & b.viewkeys()
        every_key = a.viewkeys() | b.viewkeys()
        return {k: merge(a[k], b[k]) if k in a_and_b else 
                deepcopy(a[k] if k in a else b[k]) for k in every_key}
    return deepcopy(b)


#-----------------------------------------------------------------------
def gather(filename):

    # Get local parameters
    import json
    fin = open(filename, 'r')
    dataLocal = json.load(fin)
    fin.close()

    if "include" in dataLocal.keys():

        data = {}
        for filename in dataLocal['include']:
            fin = open(filename, "r")
            dataNew = json.load(fin)
            fin.close()
            data = merge(data, dataNew)
        data = merge(data, dataLocal)
    else:
        data = dataLocal
    return data


# End of file
