# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================

import obspy


#-----------------------------------------------------------------------
def groupByStation(inventory, stream):

    stations = {}
    for network in inventory.networks:
        for station in network.stations:
            key = "%s.%s" % (network.code, station.code)
            stations[key] = []

    for trace in stream.traces:
        key = "%s.%s" % (trace.stats.network, trace.stats.station)
        stations[key].append(trace)

    for key,traces in stations.items():
        if len(traces) == 0:
            stations.pop(key)
        elif len(traces) != 3:
            print "Skipping station '%s' has %d traces. Expected 3." % (key, len(traces))
            stations.pop(key)
        else:
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
def component(stream, c):
    if c == "E":
        stA = stream.select(component="E")
        stB = stream.select(component="1")
    elif c == "N":
        stA = stream.select(component="N")
        stB = stream.select(component="2")
    elif c == "Z":
        stA = stream.select(component="Z")
        stB = stream.select(component="3")
    elif c == "R":
        return stream.select(component="R")
    elif c == "T":
        return stream.select(component="T")
    return stA + stB


# End of file
