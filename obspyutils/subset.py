# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================

import obspy


#-----------------------------------------------------------------------
def streamByStation(stream):

    stations = {}
    for trace in stream.traces:
        key = "%s.%s" % (trace.stats.network, trace.stats.station)
        if not key in stations.keys():
            stations[key] = []
        stations[key].append(trace)

    for key,traces in stations.items():
        s = obspy.core.Stream(traces=traces)
        stations[key] = s

    return stations


#-----------------------------------------------------------------------
def azimuth(stream, minimum=0.0, maximum=360.0):
    traces = []

    if maximum > minimum:
        for trace in stream:
            if trace.stats.azimuth >= minimum and trace.stats.azimuth <= maximum:
                traces.append(trace)
    else:
        for trace in stream:
            if trace.stats.azimuth >= minimum or trace.stats.azimuth <= maximum:
                traces.append(trace)

    return obspy.core.Stream(traces=traces)


#-----------------------------------------------------------------------
def pairsObsSyn(obs, syn):
    """
    Create new stream that pairs channels for observations and synthetics.
    """

    if len(obs) != len(syn):
        print "WARNING: Number of observed traces (%d) does not match number of synthetic traces (%d)." % \
                         (len(obs), len(syn))
    
    pairs = []
    for trObs in obs:
        trSyn = syn.select(network=trObs.stats.network, station=trObs.stats.station, channel="??"+trObs.stats.channel[-1]).traces
        if 1 == len(trSyn):
            pairs.append((trObs, trSyn[0]))
        else:
            print "WARNING: No synthetic found matching trace: %s" % trObs
    return pairs




# End of file
